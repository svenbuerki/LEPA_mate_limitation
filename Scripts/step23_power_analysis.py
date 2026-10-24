"""step23_power_analysis.py — can our sample detect mate limitation?

The test: do mothers with fewer compatible mates set less seed than their
plant size predicts? Seed–mother genotyping measures each mother's
compatibility; Step 13 measures her seed yield relative to her size.
This script estimates, by simulation, how often a sampling design
detects the link, for a range of effect sizes.

Simulation (per replicate):

1. **Demes.** Every deme of each population at its real size, simulated
   plant by plant (`srk_detection`: tetraploid, empirical zygosity,
   species-wide SRK frequencies).
2. **True compatibility** of a plant = share of the other plants of its
   deme it can mate with under sporophytic Class I / II SI. Plants with no
   compatible mate set no seed and are never collected.
3. **Seed yield** relative to the size expectation (log10) = effect +
   population term + plant residual. Effect: a mother with half the
   typical compatibility sets X % less seed (seeds ∝ compatibility^γ,
   γ = log(1 − X) / log 0.5; X = 10 – 50 %). Plant residual SD = Step 13's
   size → seed residual SD (0.41 log10); population term SD = Step 13's
   between-population spread beyond plant noise (computed from
   `step30i_size_seed_population_strata.tsv`). Both already contain any
   real mate-limitation signal, so the noise is if anything overstated.
4. **Sampling** (design). Mothers drawn at random among the plants that
   set seed in each sampled deme.
5. **Estimated compatibility** of a mother = share of the other genotyped
   plants of her deme she can mate with (her genotype is known from her
   seedlings); with fewer than `MIN_REFERENCE` of them, the genotyped
   plants of her whole population are used.
6. **Tests** (one-sided, α = 0.05):
   - *within populations* — mother-level regression of seed yield on
     estimated compatibility with one intercept per population and season;
   - *between populations* — population-season mean seed yield on mean
     estimated compatibility (only where there are many populations).

Scenarios:
- D1 `plan` — the Step 28 plan; `plan_all_adults` — the plan plus leaf
  genotypes of every adult in the selected demes; `all_fit_demes` — every
  fit deme of P34 / P36 at its Step 1 target.
- `all_populations` — every population sampled well enough (Step 21),
  with the mothers we hold, selected as Step 21 does (≤ Step 1 target per
  deme); also random subsets of `SUBSAMPLE` mothers from it, to show how
  many mothers the test needs.

Outputs: tables/step23_power.tsv, figures/step23_power.{png,pdf}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from srk_detection import _bits, _genotypes, _model

TABLES = Path("tables")
FIGURES = Path("figures")
DEMES_TSV = TABLES / "step30_D_model_design_demes.tsv"
STRATA_TSV = TABLES / "step30i_size_seed_population_strata.tsv"
YEARS = (2025, 2026)
RESIDUAL_SD = 0.413            # Step 13 size → seed residual SD (log10), every plant
LESS_SEED = (0.10, 0.20, 0.30, 0.50)
MIN_REFERENCE = 3
ALPHA = 0.05
N_REPS = 300
SUBSAMPLE = (100, 200, 400, 800)   # mothers drawn from the all-population sample
RNG_SEED = 2036
LABELS = {"plan": "D1 plan (selected demes)",
          "plan_all_adults": "D1 plan + every adult of the selected demes genotyped",
          "all_fit_demes": "D1, every fit deme sampled",
          "all_populations": "Every population sampled well enough"}


def between_population_sd() -> float:
    """SD of population-season mean seed yield beyond plant noise."""
    t = pd.read_csv(STRATA_TSV, sep="\t")
    t = t[t["n_plants"] >= 3]
    se2 = (t["ci95"] / 1.96) ** 2
    return float(np.sqrt(max(t["mean_resid"].var() - se2.mean(), 0.0)))


def design_tables() -> dict[str, dict]:
    """scenario → {(pop, year): DataFrame[demeID, plants, mothers]} of the
    demes it samples (mothers = 0 for demes simulated but not sampled)."""
    d1 = pd.read_csv(DEMES_TSV, sep="\t")
    out = {}
    for sc in ("plan", "plan_all_adults", "all_fit_demes"):
        units = {}
        for key, g in d1.groupby(["populationID", "year"]):
            m = (g["mothers_used"].where(g["selected"], 0) if sc != "all_fit_demes"
                 else np.minimum(g["mother_target"], g["mothers_held"]).where(g["fit"], 0))
            units[key] = pd.DataFrame(dict(demeID=g["demeID"], plants=g["fertile_plants"],
                                           mothers=m.astype(int)))
        out[sc] = units
    units = {}
    for yr in YEARS:
        pops = pd.read_csv(TABLES / f"step30_A_sampling_adequacy_{yr}.tsv", sep="\t")
        keep = set(pops.loc[pops["include_in_experiment"], "populationID"])
        dm = pd.read_csv(TABLES / f"step30_A_sampling_adequacy_per_deme_{yr}.tsv", sep="\t")
        for pop, g in dm[dm["populationID"].isin(keep)].groupby("populationID"):
            units[(pop, yr)] = pd.DataFrame(dict(demeID=g["demeID"], plants=g["deme_N_fertile"],
                                                 mothers=g["M_selected_capacity"].astype(int)))
    out["all_populations"] = units
    return out


def _expr(g: np.ndarray) -> np.ndarray:
    _, class_i, _ = _model()
    ci = np.uint64(_bits(np.where(class_i)[0]))
    own = _bits(g)
    return np.where(own & ci, own & ci, own)


def simulate_and_sample(unit: pd.DataFrame, all_adults: bool, rng):
    """Simulate one population-season and sample it. Returns arrays of the
    sampled mothers' true and estimated compatibility."""
    demes = []
    for _, r in unit.iterrows():
        n = int(r["plants"])
        if n <= 0:
            continue
        e = _expr(_genotypes(1, n, rng)[0])
        comp = (e[:, None] & e[None, :]) == 0
        np.fill_diagonal(comp, False)
        c = comp.sum(axis=1) / max(n - 1, 1)
        m = int(r["mothers"])
        seeders = np.nonzero(c > 0)[0]
        mom = (rng.choice(seeders, size=min(m, len(seeders)), replace=False)
               if m > 0 and len(seeders) else np.array([], int))
        ref = np.arange(n) if (all_adults and len(mom)) else mom
        demes.append((e, c, mom, ref))
    pop_ref = np.concatenate([e[ref] for e, _, _, ref in demes]) if demes else np.array([], np.uint64)
    ct, ce = [], []
    for e, c, mom, ref in demes:
        for i in mom:
            own = ref[ref != i]
            if len(own) >= MIN_REFERENCE:
                est = float(((e[i] & e[own]) == 0).mean())
            else:
                hit = np.nonzero(pop_ref == e[i])[0]
                pr = np.delete(pop_ref, hit[0]) if len(hit) else pop_ref
                est = float(((e[i] & pr) == 0).mean()) if len(pr) else np.nan
            ct.append(c[i]); ce.append(est)
    return np.array(ct), np.array(ce)


def within_test(y, x, grp) -> float:
    ok = np.isfinite(x)
    y, x, grp = y[ok], x[ok], grp[ok]
    X = np.column_stack([x, pd.get_dummies(grp).to_numpy(float)])
    dof = len(y) - X.shape[1]
    if dof < 2:
        return 1.0
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    cov = (res @ res / dof) * np.linalg.pinv(X.T @ X)
    if cov[0, 0] <= 0:
        return 1.0
    return float(stats.t.sf(beta[0] / np.sqrt(cov[0, 0]), dof))


def between_test(ym, xm) -> float:
    ok = np.isfinite(xm)
    if ok.sum() < 5:
        return 1.0
    r = stats.linregress(xm[ok], ym[ok])
    return float(r.pvalue / 2 if r.slope > 0 else 1 - r.pvalue / 2)


def main() -> None:
    tau = between_population_sd()
    designs = design_tables()
    gammas = {x: np.log(1 - x) / np.log(0.5) for x in LESS_SEED}
    rng = np.random.default_rng(RNG_SEED)
    scopes = {"plan": [None, 34, 36], "plan_all_adults": [None], "all_fit_demes": [None],
              "all_populations": [None]}
    hits, n_mothers, n_units, spread = {}, {}, {}, []
    for sc, units in designs.items():
        for rep in range(N_REPS):
            sampled = {k: simulate_and_sample(u, sc == "plan_all_adults", rng)
                       for k, u in units.items()}
            sampled = {k: v for k, v in sampled.items() if len(v[0])}
            if rep == 0:
                n_units[sc] = len(sampled)
                for s in scopes[sc]:
                    n_mothers[(sc, s)] = sum(len(v[0]) for k, v in sampled.items()
                                             if s is None or k[0] == s)
            for k, (ct, _) in sampled.items():
                if rep < 20:
                    spread.append(dict(scenario=sc, pop=k[0], year=k[1], mean=ct.mean(), sd=ct.std()))
            plant_noise = {k: rng.normal(0, RESIDUAL_SD, len(v[0])) for k, v in sampled.items()}
            pop_term = {k: rng.normal(0, tau) for k in sampled}
            for x, gm in gammas.items():
                # one common reference compatibility (species-wide mean of seed-setters)
                for s in scopes[sc]:
                    ys, xs, gs, ym, xm = [], [], [], [], []
                    for k, (ct, ce) in sampled.items():
                        if s is not None and k[0] != s:
                            continue
                        y = gm * np.log10(ct / 0.78) + pop_term[k] + plant_noise[k]
                        ys.append(y); xs.append(ce); gs.append(np.full(len(ct), f"{k}"))
                        ym.append(y.mean())
                        xm.append(np.nanmean(ce) if np.isfinite(ce).any() else np.nan)
                    if not ys:
                        continue
                    p_w = within_test(np.concatenate(ys), np.concatenate(xs), np.concatenate(gs))
                    hits[(sc, s, "within", x)] = hits.get((sc, s, "within", x), 0) + (p_w < ALPHA)
                    if sc == "all_populations":
                        p_b = between_test(np.array(ym), np.array(xm))
                        hits[(sc, s, "between", x)] = hits.get((sc, s, "between", x), 0) + (p_b < ALPHA)
                        Y, Xc, Gc = (np.concatenate(ys), np.concatenate(xs), np.concatenate(gs))
                        for K in SUBSAMPLE:
                            if K >= len(Y):
                                continue
                            idx = rng.choice(len(Y), size=K, replace=False)
                            p_k = within_test(Y[idx], Xc[idx], Gc[idx])
                            hits[(sc, K, "within", x)] = hits.get((sc, K, "within", x), 0) + (p_k < ALPHA)
                            n_mothers[(sc, K)] = K
        print(f"[step23] {sc}: done")
    rows = []
    for (sc, s, test, x), h in hits.items():
        scope = ("all" if s is None else f"P{s}" if sc != "all_populations"
                 else f"subset of {s}")
        rows.append(dict(scenario=sc, label=LABELS[sc], scope=scope, test=test,
                         less_seed_at_half_compatibility=x, mothers=n_mothers[(sc, s)],
                         population_seasons=n_units[sc] if s is None or sc == "all_populations" else None,
                         power=round(h / N_REPS, 3), reps=N_REPS))
    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "step23_power.tsv", sep="\t", index=False)
    sp = pd.DataFrame(spread).groupby("scenario")[["mean", "sd"]].mean().round(3)
    print(f"[step23] between-population SD of seed yield beyond plant noise: {tau:.3f} log10")
    print("[step23] true compatibility of sampled mothers (mean, SD within a population):\n",
          sp.to_string())
    print(out.pivot_table(index=["scenario", "scope", "test"],
                          columns="less_seed_at_half_compatibility", values="power").to_string())
    plot(out)


def plot(out: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4))
    ax = axes[0]
    style = {("plan", "all"): ("#E41A1C", "-"), ("plan", "P34"): ("#FF7F00", "--"),
             ("plan", "P36"): ("#FDBF6F", "--"), ("plan_all_adults", "all"): ("#984EA3", "-"),
             ("all_fit_demes", "all"): ("#377EB8", "-")}
    d1 = out[out["scenario"] != "all_populations"]
    for (sc, scope), g in d1.groupby(["scenario", "scope"]):
        c, ls = style[(sc, scope)]
        g = g.sort_values("less_seed_at_half_compatibility")
        name = LABELS[sc] if scope == "all" else f"D1 plan, {scope} only"
        ax.plot(100 * g["less_seed_at_half_compatibility"], g["power"], "o", ls=ls, color=c,
                label=f"{name} — {int(g['mothers'].iat[0])} mothers")
    ax.set_xlabel("How much less seed a mother sets with half the usual share of compatible mates (%)")
    ax.set_ylabel("Chance of detecting mate limitation")
    ax.set_title("A. D1 model populations — designs compared", loc="left", weight="bold")
    ax = axes[1]
    ap = out[(out["scenario"] == "all_populations") & (out["test"] == "within")].copy()
    ap["n"] = ap["mothers"]
    cmap = plt.get_cmap("viridis")
    for k, (x, g) in enumerate(ap.groupby("less_seed_at_half_compatibility")):
        g = g.sort_values("n")
        ax.plot(g["n"], g["power"], "o-", color=cmap(k / 3), label=f"{100 * x:.0f} % less seed")
    allp = ap[ap["scope"] == "all"]["n"].iat[0]
    ax.axvline(allp, ls="--", c="#888", lw=0.8)
    ax.text(allp, 0.03, f" all {allp} mothers held\n ({int(ap['population_seasons'].iat[0])} population-seasons)",
            fontsize=7.5, color="#555")
    ax.axvline(91, ls="--", c="#E41A1C", lw=0.8)
    ax.text(91, 0.03, " D1 plan size", fontsize=7.5, color="#E41A1C")
    ax.set_xscale("log")
    ax.set_xlabel("Mothers genotyped, drawn across every population sampled well enough")
    ax.set_title("B. How many mothers does the test need?", loc="left", weight="bold")
    ax.legend(title="Effect: with half the compatible mates", fontsize=7.5, title_fontsize=7.5,
              loc="upper left")
    for ax in axes:
        ax.axhline(0.8, ls=":", c="#444", lw=0.8)
        ax.set_ylim(0, 1.02)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(fontsize=7.5, loc="upper left")
    fig.suptitle("Can our sample detect mate limitation? (simulation; dotted line = 80 %)",
                 weight="bold")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"step23_power.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
