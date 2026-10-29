"""step33_srk_drift.py — SRK-specific drift analysis at the genotyped sites
(part 1: allele number, evenness, sharing).

The calibration of Step 31 treats SRK as a neutral gene. SRK is not:
self-incompatibility favours rare alleles, and in sporophytic SI with a
dominance hierarchy recessive alleles (Class II in LEPA) reach higher
frequencies than dominant ones (Schierup, Vekemans & Christiansen 1997,
Genetics 147: 835 – 846). This step asks three questions of the six
genotyped sites (the plants of Step 25; genotyped plants treated as
independent):

1. **How many SRK alleles does each population hold?** Samples are small
   (12 – 77 plants), so some alleles are missed. The incidence-based Chao2
   estimator (Chao 1987, Biometrics 43: 783 – 791) adds the alleles likely
   unseen, from those found in exactly one and two plants; 95 % interval
   from Chao's log-normal formula (it never falls below the alleles seen).
2. **Are the alleles more unequal than the SI system predicts?** Under
   sporophytic dominance, alleles of the same class are expected at about
   equal frequency. For each class, the observed unevenness (Simpson's
   concentration of the class's copies) is compared with plants simulated
   from equal within-class frequencies (observed class shares and number of
   alleles, empirical zygosity, same number of plants). Excess unevenness
   beyond sampling is a drift signal; Class II being commoner than Class I
   is not.
3. **Which alleles do the sites share?** Jaccard similarity of allele sets.

Part 2 — forward simulation of tetraploid sporophytic SI. A population of
N tetraploid plants is founded from the species-wide alleles (empirical
zygosity). Each generation, every offspring has a random mother and a
father drawn among the plants compatible with her (Class I / II
dominance, expressed sets disjoint); each parent passes 2 of its 4 copies
(random chromosome segregation). After T generations, as many plants as
were genotyped are sampled. The same run with any father allowed is the
neutral comparison. For each site and T, the N reproducing the observed
number of SRK alleles best (likelihood of the exact count over N_SIM_FWD
runs) is its effective size under SI and under neutrality.

Outputs: tables/step33_srk_alleles_per_site.tsv,
tables/step33_si_simulation_grid.tsv, tables/step33_si_simulation_fit.tsv,
tables/step33_srk_evenness.tsv, tables/step33_srk_sharing.tsv,
figures/step33_srk_drift.{png,pdf}, figures/step33_si_simulation.{png,pdf}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TABLES = Path("tables")
FIGURES = Path("figures")
SITES = {"EO67": "P6 (EO67)", "EO70": "P42 (EO70)", "EO76": "P40 + P41 (EO76)",
         "P34": "P34 (EO18)", "P13": "P13 (EO27)", "P16": "P16 (EO27)"}
N_BOOT = 1000
N_NULL = 2000
NE_GRID = (5, 8, 12, 18, 25, 35, 50, 75, 100, 150, 250, 400)
T_GRID = (10, 25, 50)
N_SIM_FWD = 120
RNG_SEED = 2039
plt.rcParams.update({"font.size": 13})


def load_sites() -> tuple[dict, list[str], np.ndarray]:
    """Per site: (plants × alleles) copy-count matrix of functional SRK copies."""
    import step25_partC_clean_overlap as s25
    from step24_EO_scale_validation import load_allele_to_fg, load_individuals
    from srk_si_model import build_class_i_mask, load_class_map
    prior = pd.read_csv(s25.DEFAULT_PRIOR_TSV, sep="\t", encoding="utf-8-sig")
    fg = prior["Fg"].astype(str).tolist()
    ind = load_individuals(s25.DEFAULT_INDIV_TSV, load_allele_to_fg(s25.DEFAULT_ALLELE_FG_TSV), fg)
    pc = pd.read_csv(s25.DEFAULT_PRED_PCOMPAT, sep="\t", encoding="utf-8-sig")
    dv = pd.read_csv(s25.DEFAULT_PRED_DIVERSITY, sep="\t", encoding="utf-8-sig")
    ind, _, _, _ = s25.add_population_sites(ind, pc, dv)
    out = {}
    for site in SITES:
        sub = ind[(ind["EO"] == site) & (ind["n_functional_copies"] > 0)]
        out[site] = np.stack(sub["fg_copies"].values)
    return out, fg, build_class_i_mask(fg, load_class_map())


def chao2(inc: np.ndarray) -> float:
    """Incidence-based richness estimate from a (plants × alleles) 0/1 matrix."""
    m = inc.shape[0]
    q = inc.sum(axis=0)
    s_obs = int((q > 0).sum())
    f1, f2 = int((q == 1).sum()), int((q == 2).sum())
    if f2 > 0:
        return s_obs + (m - 1) / m * f1 ** 2 / (2 * f2)
    return s_obs + (m - 1) / m * f1 * (f1 - 1) / 2


def chao2_ci(inc: np.ndarray) -> tuple[float, float]:
    """Chao's (1987) log-normal 95 % interval for the Chao2 estimate."""
    m = inc.shape[0]
    q = inc.sum(axis=0)
    s_obs = int((q > 0).sum())
    f1, f2 = int((q == 1).sum()), int((q == 2).sum())
    a = (m - 1) / m
    if f2 > 0:
        r = f1 / f2
        var = f2 * (a * r ** 2 / 2 + a ** 2 * r ** 3 + a ** 2 * r ** 4 / 4)
    else:
        var = (a * f1 * (f1 - 1) / 2 + a ** 2 * f1 * (2 * f1 - 1) ** 2 / 4
               - a ** 2 * f1 ** 4 / (4 * max(chao2(inc), 1e-9)))
    t = chao2(inc) - s_obs
    if t <= 0 or var <= 0:
        return float(s_obs), float(s_obs)
    k = np.exp(1.96 * np.sqrt(np.log(1 + var / t ** 2)))
    return s_obs + t / k, s_obs + t * k


def richness(sites: dict, rng) -> pd.DataFrame:
    rows = []
    for site, x in sites.items():
        inc = (x > 0).astype(int)
        lo, hi = chao2_ci(inc)
        q = inc.sum(axis=0)
        rows.append(dict(site=site, label=SITES[site], plants=len(inc),
                         alleles_observed=int((q > 0).sum()),
                         in_one_plant=int((q == 1).sum()), in_two_plants=int((q == 2).sum()),
                         alleles_estimated=round(chao2(inc), 1),
                         estimate_lo95=round(lo, 1), estimate_hi95=round(hi, 1)))
    return pd.DataFrame(rows)


def simpson(c: np.ndarray) -> float:
    """Concentration of copies among alleles: Σ p², 1/k when all equal."""
    s = c.sum()
    return float(((c / s) ** 2).sum()) if s > 0 else np.nan


def evenness(sites: dict, class_i: np.ndarray, rng) -> pd.DataFrame:
    """Within-class unevenness against equal within-class frequencies."""
    from srk_detection import _genotypes
    rows = []
    for site, x in sites.items():
        tot = x.sum(axis=0).astype(float)
        f_obs = tot / tot.sum()
        for cls, mask in (("Class I", class_i), ("Class II", ~class_i)):
            present = (tot > 0) & mask
            k = int(present.sum())
            obs = simpson(tot[present]) if k > 1 else np.nan
            share = float(f_obs[mask].sum())
            if k < 2:
                rows.append(dict(site=site, label=SITES[site], allele_class=cls, alleles=k,
                                 class_share=round(share, 3), concentration_observed=obs,
                                 concentration_null_mean=np.nan, p_more_uneven=np.nan))
                continue
            # null: same class shares and alleles, equal frequencies within each class
            f0 = np.zeros(len(tot))
            for cm in (class_i, ~class_i):
                pres = (tot > 0) & cm
                if pres.any():
                    f0[pres] = f_obs[cm].sum() / pres.sum()
            null = np.empty(N_NULL)
            for b in range(N_NULL):
                g = _genotypes(1, len(x), rng, f0)[0]
                cnt = np.bincount(g.ravel(), minlength=len(tot)).astype(float)
                null[b] = simpson(cnt[present])
            rows.append(dict(site=site, label=SITES[site], allele_class=cls, alleles=k,
                             class_share=round(share, 3), concentration_observed=round(obs, 3),
                             concentration_null_mean=round(float(np.nanmean(null)), 3),
                             p_more_uneven=round(float(np.mean(null >= obs)), 4)))
    return pd.DataFrame(rows)


def sharing(sites: dict) -> pd.DataFrame:
    sets = {s: set(np.nonzero(x.sum(axis=0))[0]) for s, x in sites.items()}
    names = list(sets)
    mat = pd.DataFrame(index=[SITES[s] for s in names], columns=[SITES[s] for s in names], dtype=float)
    for a in names:
        for b in names:
            mat.loc[SITES[a], SITES[b]] = len(sets[a] & sets[b]) / len(sets[a] | sets[b])
    return mat


def plot(rich: pd.DataFrame, even: pd.DataFrame, share: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(19, 6.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.15, 1.05], wspace=0.45)
    ax = fig.add_subplot(gs[0])
    y = np.arange(len(rich))[::-1]
    ax.errorbar(rich["alleles_estimated"], y,
                xerr=[rich["alleles_estimated"] - rich["estimate_lo95"],
                      rich["estimate_hi95"] - rich["alleles_estimated"]],
                fmt="s", color="#377EB8", capsize=4, label="estimated in the population (95 %)")
    ax.scatter(rich["alleles_observed"], y, color="black", zorder=3, s=50, label="found in the sample")
    ax.set_yticks(y); ax.set_yticklabels([f"{l}\n{p} plants" for l, p in zip(rich["label"], rich["plants"])])
    ax.set_xlim(0, 32); ax.set_xlabel("SRK alleles (species: 32)")
    ax.set_ylim(-0.8, len(rich) - 0.2)
    ax.legend(fontsize=10, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=1, frameon=False)
    ax.set_title("A. How many SRK alleles?", loc="left", weight="bold")
    ax = fig.add_subplot(gs[1])
    labs = list(SITES.values())
    for k, cls in enumerate(("Class I", "Class II")):
        e = even[even["allele_class"] == cls].set_index("label").reindex(labs)
        yy = np.arange(len(labs))[::-1] + (0.15 if k == 0 else -0.15)
        col = "#984EA3" if k == 0 else "#4D4D4D"
        ax.scatter(e["concentration_null_mean"], yy, marker="|", s=300, color=col, alpha=0.6)
        sig = e["p_more_uneven"] < 0.05
        ax.scatter(e["concentration_observed"], yy, s=70, color=col,
                   facecolors=np.where(sig, col, "white"), edgecolors=col, linewidths=2,
                   label=f"{cls} alleles", zorder=3)
    ax.set_yticks(np.arange(len(labs))[::-1]); ax.set_yticklabels(labs)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Concentration within the class\n(| = expected if its alleles were equally common)")
    ax.legend(fontsize=10, loc="lower right", title="filled = more uneven\nthan sampling (p < 0.05)",
              title_fontsize=9)
    ax.set_title("B. Are alleles more unequal than SI predicts?", loc="left", weight="bold")
    ax = fig.add_subplot(gs[2])
    im = ax.imshow(share.values.astype(float), cmap="Greys", vmin=0, vmax=1)
    ax.set_xticks(range(len(share))); ax.set_xticklabels(share.columns, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(len(share))); ax.set_yticklabels(share.index, fontsize=10)
    for i in range(len(share)):
        for j in range(len(share)):
            v = share.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=9,
                    color="white" if v > 0.55 else "black")
    fig.colorbar(im, ax=ax, fraction=0.046, label="Shared alleles (Jaccard)")
    ax.set_title("C. Which alleles do sites share?", loc="left", weight="bold")
    for a in fig.axes[:2]:
        a.spines[["top", "right"]].set_visible(False)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"step33_srk_drift.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- part 2
def _expressed(own, ci_bits):
    return np.where(own & ci_bits, own & ci_bits, own)


def forward(N: int, T: int, n_sample: int, si: bool, rng) -> tuple[int, float]:
    """One forward run; returns (alleles in a sample of n plants, Class I
    concentration in that sample)."""
    from srk_detection import _bits, _genotypes, _model
    _, class_i, _ = _model()
    ci_bits = np.uint64(_bits(np.where(class_i)[0]))
    g = _genotypes(1, N, rng)[0]                                  # (N, 4)
    pairs = np.array([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
    for _ in range(T):
        own = _bits(g)
        ex = _expressed(own, ci_bits)
        mom = rng.integers(0, N, N)
        if si:
            dad = np.full(N, -1)
            todo = np.arange(N)
            for _try in range(60):
                if not len(todo):
                    break
                cand = rng.integers(0, N, len(todo))
                ok = (ex[cand] & ex[mom[todo]]) == 0
                dad[todo[ok]] = cand[ok]
                todo = todo[~ok]
                if len(todo):                                     # re-draw mothers left without a mate
                    mom[todo] = rng.integers(0, N, len(todo))
            if len(todo):                                         # no compatible pair at all: population fails
                return 0, np.nan
        else:
            dad = rng.integers(0, N, N)
        pm = pairs[rng.integers(0, 6, N)]
        pd_ = pairs[rng.integers(0, 6, N)]
        r = np.arange(N)[:, None]
        g = np.concatenate([g[mom][r, pm], g[dad][r, pd_]], axis=1)
    s = g[rng.integers(0, N, n_sample)]
    cnt = np.bincount(s.ravel(), minlength=len(class_i)).astype(float)
    ci = cnt[class_i]
    conc = simpson(ci[ci > 0]) if (ci > 0).sum() > 1 else np.nan
    return int((cnt > 0).sum()), conc


def si_simulation(rich: pd.DataFrame, rng) -> tuple[pd.DataFrame, pd.DataFrame]:
    sizes = sorted(set(rich["plants"]))
    rows = []
    for si in (True, False):
        for T in T_GRID:
            for N in NE_GRID:
                for n in sizes:
                    res = [forward(N, T, n, si, rng) for _ in range(N_SIM_FWD)]
                    a = np.array([x[0] for x in res])
                    rows.append(dict(model="SI" if si else "neutral", generations=T, Ne=N,
                                     sample_plants=n, alleles_mean=a.mean(),
                                     alleles_lo=np.quantile(a, .025), alleles_hi=np.quantile(a, .975),
                                     failed_share=float(np.mean(a == 0)),
                                     counts=",".join(map(str, np.bincount(a, minlength=33)))))
        print(f"[step33] simulated {'SI' if si else 'neutral'}")
    grid = pd.DataFrame(rows)
    fits = []
    for _, s in rich.iterrows():
        for (model, T), g in grid[grid["sample_plants"] == s["plants"]].groupby(["model", "generations"]):
            lik = []
            for _, r in g.iterrows():
                c = np.array(r["counts"].split(","), int)
                lik.append((c[s["alleles_observed"]] + 0.5) / (c.sum() + 16.5))
            lik = np.array(lik)
            best = g["Ne"].to_numpy()[np.argmax(lik)]
            ok = g["Ne"].to_numpy()[lik >= lik.max() / 8]            # ~95 % support (log-lik within 2)
            fits.append(dict(site=s["site"], label=s["label"], model=model, generations=T,
                             alleles_observed=s["alleles_observed"], Ne_best=int(best),
                             Ne_lo=int(ok.min()), Ne_hi=int(ok.max())))
    return grid.drop(columns="counts"), pd.DataFrame(fits)


def plot_simulation(grid: pd.DataFrame, fits: pd.DataFrame, rich: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.2))
    ax = axes[0]
    n_ref = int(rich["plants"].median())
    n_ref = min(set(rich["plants"]), key=lambda v: abs(v - n_ref))
    for model, col in (("SI", "#111111"), ("neutral", "#9E9E9E")):
        g = grid[(grid["model"] == model) & (grid["generations"] == 25) & (grid["sample_plants"] == n_ref)]
        ax.plot(g["Ne"], g["alleles_mean"], "o-", color=col, lw=2,
                label="with self-incompatibility" if model == "SI" else "neutral (no SI)")
        ax.fill_between(g["Ne"], g["alleles_lo"], g["alleles_hi"], color=col, alpha=0.15)
    ax.axhspan(rich["alleles_observed"].min(), rich["alleles_observed"].max(), color="#8FB8DE",
               alpha=0.15, label="observed at the six sites")
    ax.set_xscale("log"); ax.set_xlabel("Effective population size (plants)")
    ax.set_ylabel(f"SRK alleles in a sample of {n_ref} plants")
    ax.set_title("A. Alleles kept after 25 generations", loc="left", weight="bold")
    ax.legend(fontsize=11, loc="lower right")
    ax = axes[1]
    f25 = fits[fits["generations"] == 25]
    labs = list(SITES.values())
    for k, (model, col) in enumerate((("SI", "#111111"), ("neutral", "#9E9E9E"))):
        d = f25[f25["model"] == model].set_index("label").reindex(labs)
        yy = np.arange(len(labs))[::-1] + (0.14 if k == 0 else -0.14)
        ax.errorbar(d["Ne_best"], yy, xerr=[d["Ne_best"] - d["Ne_lo"], d["Ne_hi"] - d["Ne_best"]],
                    fmt="o", color=col, capsize=4, label="with self-incompatibility" if model == "SI" else "neutral (no SI)")
    ax.set_yticks(np.arange(len(labs))[::-1]); ax.set_yticklabels(labs)
    ax.set_xscale("log"); ax.set_xlabel("Effective population size implied (plants)")
    ax.set_title("B. Effective size each site implies (25 generations)", loc="left", weight="bold")
    ax.set_ylim(-0.7, len(labs) - 0.3)
    ax.legend(fontsize=11, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=False)
    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"step33_si_simulation.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)
    sites, fg, class_i = load_sites()
    rich = richness(sites, rng)
    even = evenness(sites, class_i, rng)
    share = sharing(sites)
    rich.to_csv(TABLES / "step33_srk_alleles_per_site.tsv", sep="\t", index=False)
    even.to_csv(TABLES / "step33_srk_evenness.tsv", sep="\t", index=False)
    share.round(3).to_csv(TABLES / "step33_srk_sharing.tsv", sep="\t")
    pd.set_option("display.width", 220)
    print(rich.to_string(index=False)); print(even.to_string(index=False)); print(share.round(2).to_string())
    plot(rich, even, share)
    grid, fits = si_simulation(rich, rng)
    grid.to_csv(TABLES / "step33_si_simulation_grid.tsv", sep="\t", index=False)
    fits.to_csv(TABLES / "step33_si_simulation_fit.tsv", sep="\t", index=False)
    print(fits.pivot_table(index="label", columns=["model", "generations"], values="Ne_best").to_string())
    plot_simulation(grid, fits, rich)
    print("[step33] wrote tables/step33_*.tsv and figures/step33_*.{png,pdf}")


if __name__ == "__main__":
    main()
