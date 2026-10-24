"""step31_drift_calibration.py — calibrate the drift model on the
SRK-genotyped sites.

The drift model of the predictions draws each deme from the species-wide
SRK frequencies. At five of the six genotyped sites it predicts 22 – 29
alleles where 5 – 9 are observed (Step 25): real populations have drifted
much further. This step fits how far.

Model. A population's SRK pool descends from **Ne effective founder
plants**: Ne tetraploid plants (empirical zygosity) drawn from the
species-wide frequencies set the population's allele frequencies; its
plants — and the plants we genotyped — are drawn from that pool. Ne is
tied to the census N by

    Ne = a · N^b        (rounded, ≥ 1)

so b = 1 is "Ne a fixed share of the census" and b = 0 is "every
population descends from the same number of founders, whatever its size
today". Genotyped plants are treated as independent (no family
correction; Sven, 2026-10-24).

Fit. For each (a, b) on a grid and each site, the number of distinct SRK
alleles among n genotyped plants is simulated (`N_SIM` populations); the
likelihood is the simulated probability of the observed count. The best
(a, b) maximises the sum over sites. Checks:
- **Leave one site out** — fit on five sites, predict the sixth.
- **Pollen compatibility** (not used in the fit) — simulated with the same
  statistic as Step 25 (each sampled plant against fathers drawn from the
  sample's allele frequencies with the empirical zygosity).

Re-prediction. With the fitted (a, b), every population-season is
re-predicted: SRK alleles among its N plants and mean pollen
compatibility, and compared with the seed-yield read-out (Step 13). The
primary predictions (Step 02) are left unchanged.

Inputs: Phase 4 genotypes + placement (as Step 25) ·
tables/step30_B_partC_clean_overlap_per_location.tsv ·
tables/step30g_prediction_population_year.tsv ·
tables/step30i_seed_yield_vs_predictions.tsv
Outputs: tables/step31_drift_calibration_fit.tsv,
tables/step31_drift_calibration_sites.tsv,
tables/step31_calibrated_predictions.tsv,
figures/step31_drift_calibration.{png,pdf}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from srk_bl_constants import BL_COLORS
from srk_detection import _bits, _genotypes, _model, deme_allele_counts, popcount
from srk_si_model import load_zygosity_dist, p_compat_sporophytic_empirical

TABLES = Path("tables")
FIGURES = Path("figures")
SITES = {"EO67": [6], "EO70": [42], "EO76": [40, 41], "P34": [34], "P13": [13], "P16": [16]}
SITE_LABEL = {"EO67": "P6 (EO67)", "EO70": "P42 (EO70)", "EO76": "P40 + P41 (EO76)",
              "P34": "P34 (EO18)", "P13": "P13 (EO27)", "P16": "P16 (EO27)"}
CENSUS_YEAR = 2025                 # the genotyped seed is mostly 2025 seed
A_GRID = np.logspace(-2, 2.5, 37)
B_GRID = np.round(np.arange(0, 1.01, 0.1), 2)
N_SIM = 1000
N_SIM_COMPAT = 150
RNG_SEED = 2037


def site_table() -> pd.DataFrame:
    o = pd.read_csv(TABLES / "step30_B_partC_clean_overlap_per_location.tsv", sep="\t")
    p = pd.read_csv(TABLES / "step30g_prediction_population_year.tsv", sep="\t")
    p = p[p["year"] == CENSUS_YEAR].set_index("populationID")["N_fertile_total"]
    rows = []
    for s, pops in SITES.items():
        r = o[o["locationCode"] == s].iloc[0]
        rows.append(dict(site=s, label=SITE_LABEL[s], N=int(sum(p.get(q, 0) for q in pops)),
                         n=int(r["n_functional_carriers"]), obs_alleles=int(r["obs_distinct_Fgs"]),
                         obs_compat=float(r["obs_pcompat_mean"]),
                         obs_compat_lo=float(r["obs_pcompat_lo95"]),
                         obs_compat_hi=float(r["obs_pcompat_hi95"]),
                         uncalibrated_alleles=float(r["pred_phase5_distinct_Fgs"])))
    return pd.DataFrame(rows)


NE_CAP = 500   # beyond ~500 founders the pool is the species-wide one


def ne_of(N, a, b) -> int:
    return int(min(NE_CAP, max(1, round(a * float(N) ** b))))


_cache: dict = {}


def sample_alleles(ne: int, n: int, rng) -> np.ndarray:
    """Distinct alleles among n plants drawn from a pool founded by ne
    plants — N_SIM simulated populations (cached by (ne, n))."""
    key = (ne, n)
    if key not in _cache:
        out = np.empty(N_SIM, int)
        for r in range(N_SIM):
            f = deme_allele_counts(ne, rng)
            g = _genotypes(1, n, rng, f)[0]
            out[r] = int(popcount(np.bitwise_or.reduce(_bits(g)))[()])
        _cache[key] = out
    return _cache[key]


def loglik(sites: pd.DataFrame, a: float, b: float, rng) -> np.ndarray:
    ll = []
    for _, s in sites.iterrows():
        d = sample_alleles(ne_of(s["N"], a, b), int(s["n"]), rng)
        p = ((d == s["obs_alleles"]).sum() + 0.5) / (len(d) + 0.5 * 33)
        ll.append(np.log(p))
    return np.array(ll)


def fit(sites: pd.DataFrame, rng) -> tuple[pd.DataFrame, float, float]:
    rows = []
    for b in B_GRID:
        for a in A_GRID:
            ll = loglik(sites, a, b, rng)
            rows.append(dict(a=a, b=b, loglik=float(ll.sum())))
    g = pd.DataFrame(rows)
    best = g.loc[g["loglik"].idxmax()]
    return g, float(best["a"]), float(best["b"])


def compat_sample(ne: int, n: int, rng) -> np.ndarray:
    """Step 25's compatibility statistic on n plants from a pool founded
    by ne plants."""
    _, class_i, _ = _model()
    zyg = load_zygosity_dist()
    out = []
    for _ in range(N_SIM_COMPAT):
        f = deme_allele_counts(ne, rng)
        g = _genotypes(1, n, rng, f)[0]
        local = np.bincount(g.ravel(), minlength=len(f)) / g.size
        out.append(p_compat_sporophytic_empirical(g, local, class_i, zyg,
                                                  n_fathers=500, rng=rng).mean())
    return np.array(out)


def population_prediction(N: int, a: float, b: float, rng, reps: int = 120) -> dict:
    """Calibrated SRK alleles among the N plants and mean compatibility."""
    _, class_i, _ = _model()
    zyg = load_zygosity_dist()
    ne = ne_of(N, a, b)
    div, comp = [], []
    for _ in range(reps):
        f = deme_allele_counts(ne, rng)
        g = _genotypes(1, max(int(N), 1), rng, f)[0]
        div.append(int(popcount(np.bitwise_or.reduce(_bits(g)))[()]))
        local = np.bincount(g.ravel(), minlength=len(f)) / g.size
        moms = g[rng.choice(len(g), size=min(len(g), 200), replace=False)]
        comp.append(p_compat_sporophytic_empirical(moms, local, class_i, zyg,
                                                   n_fathers=400, rng=rng).mean())
    div, comp = np.array(div), np.array(comp)
    return dict(Ne=ne, cal_alleles_mean=div.mean(), cal_alleles_lo=np.quantile(div, .025),
                cal_alleles_hi=np.quantile(div, .975), cal_compat_mean=comp.mean(),
                cal_compat_lo=np.quantile(comp, .025), cal_compat_hi=np.quantile(comp, .975))


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)
    sites = site_table()
    grid, a, b = fit(sites, rng)
    grid.to_csv(TABLES / "step31_drift_calibration_fit.tsv", sep="\t", index=False)
    prof = grid.groupby("b")["loglik"].max()
    print(f"[step31] best fit: Ne = {a:.3g} · N^{b:.1f}; log-likelihood {grid['loglik'].max():.2f}")
    print("[step31] best log-likelihood per b:\n", prof.round(2).to_string())

    rows = []
    for i, s in sites.iterrows():
        rest = sites.drop(index=i)
        _, a_i, b_i = fit(rest, rng)
        ne_i = ne_of(s["N"], a_i, b_i)
        d_i = sample_alleles(ne_i, int(s["n"]), rng)
        ne = ne_of(s["N"], a, b)
        d = sample_alleles(ne, int(s["n"]), rng)
        c = compat_sample(ne, int(s["n"]), rng)
        rows.append(dict(**s.to_dict(), Ne=ne,
                         fit_alleles_mean=d.mean(), fit_alleles_lo=np.quantile(d, .025),
                         fit_alleles_hi=np.quantile(d, .975),
                         loo_a=a_i, loo_b=b_i, loo_Ne=ne_i, loo_alleles_mean=d_i.mean(),
                         loo_alleles_lo=np.quantile(d_i, .025), loo_alleles_hi=np.quantile(d_i, .975),
                         loo_inside=bool(np.quantile(d_i, .025) <= s["obs_alleles"] <= np.quantile(d_i, .975)),
                         fit_compat_mean=c.mean(), fit_compat_lo=np.quantile(c, .025),
                         fit_compat_hi=np.quantile(c, .975),
                         compat_inside=bool(np.quantile(c, .025) <= s["obs_compat"] <= np.quantile(c, .975))))
    st = pd.DataFrame(rows)
    st.to_csv(TABLES / "step31_drift_calibration_sites.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(st[["label", "N", "n", "Ne", "obs_alleles", "uncalibrated_alleles", "fit_alleles_mean",
              "loo_alleles_lo", "loo_alleles_hi", "loo_inside", "obs_compat", "fit_compat_mean",
              "fit_compat_lo", "fit_compat_hi", "compat_inside"]].round(2).to_string(index=False))

    pred = pd.read_csv(TABLES / "step30g_prediction_population_year.tsv", sep="\t")
    out = []
    for _, r in pred.iterrows():
        out.append(dict(populationID=int(r["populationID"]), BL=r["BL"], year=int(r["year"]),
                        N_fertile_total=int(r["N_fertile_total"]),
                        uncal_alleles_mean=r["pred_srk_diversity_mean"],
                        uncal_compat_mean=r["pred_pcompat_mean"],
                        **population_prediction(int(r["N_fertile_total"]), a, b, rng)))
    cal = pd.DataFrame(out)
    cal.to_csv(TABLES / "step31_calibrated_predictions.tsv", sep="\t", index=False)

    sy = pd.read_csv(TABLES / "step30i_seed_yield_vs_predictions.tsv", sep="\t")
    m = sy.merge(cal, on=["populationID", "year"])
    tests = []
    for col, lab in (("cal_compat_mean", "calibrated compatibility"),
                     ("cal_alleles_mean", "calibrated SRK alleles"),
                     ("uncal_compat_mean", "uncalibrated compatibility")):
        rho, p = stats.spearmanr(m[col], m["log10_fold"])
        tests.append(dict(predictor=lab, n=len(m), spearman_rho=rho, p_value=p))
    tests = pd.DataFrame(tests)
    print("[step31] seed yield (fold of expectation) vs predictions:\n", tests.round(3).to_string(index=False))
    tests.to_csv(TABLES / "step31_calibrated_vs_seed_yield.tsv", sep="\t", index=False)
    plot(grid, st, cal, m, a, b)


def plot(grid, st, cal, m, a, b) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 10.5))
    ax = axes[0, 0]
    prof = grid.groupby("b")["loglik"].max()
    ax.plot(prof.index, prof.values, "o-", c="#333")
    ax.set_xlabel("b in Ne = a · N^b   (0 = same founders whatever the size, 1 = share of census)")
    ax.set_ylabel("Best log-likelihood over a")
    ax.set_title("A. How effective size scales with census", loc="left", weight="bold")
    ax = axes[0, 1]
    y = np.arange(len(st))[::-1]
    ax.scatter(st["uncalibrated_alleles"], y, marker="x", c="#999", label="uncalibrated prediction")
    ax.errorbar(st["loo_alleles_mean"], y, xerr=[st["loo_alleles_mean"] - st["loo_alleles_lo"],
                st["loo_alleles_hi"] - st["loo_alleles_mean"]], fmt="s", c="#377EB8", capsize=3,
                label="calibrated, site left out of the fit (95 %)")
    ax.scatter(st["obs_alleles"], y, c="#E41A1C", zorder=5, label="observed")
    ax.set_yticks(y); ax.set_yticklabels(st["label"])
    ax.set_xlabel("SRK alleles among the genotyped plants")
    ax.set_title("B. SRK diversity at the genotyped sites", loc="left", weight="bold")
    ax.legend(fontsize=8, loc="lower right")
    ax = axes[1, 0]
    ax.errorbar(st["fit_compat_mean"], y, xerr=[st["fit_compat_mean"] - st["fit_compat_lo"],
                st["fit_compat_hi"] - st["fit_compat_mean"]], fmt="s", c="#377EB8", capsize=3,
                label="calibrated prediction (95 %)")
    ax.errorbar(st["obs_compat"], y, xerr=[st["obs_compat"] - st["obs_compat_lo"],
                st["obs_compat_hi"] - st["obs_compat"]], fmt="o", c="#E41A1C", capsize=3,
                label="observed (95 %)")
    ax.set_yticks(y); ax.set_yticklabels(st["label"])
    ax.set_xlabel("Pollen compatibility (share of plants a plant can mate with)")
    ax.set_title("C. Compatibility — not used in the fit", loc="left", weight="bold")
    ax.legend(fontsize=8, loc="lower left")
    ax = axes[1, 1]
    for bl, g in m.groupby("BL_x" if "BL_x" in m else "BL"):
        ax.scatter(g["cal_compat_mean"], g["log10_fold"], c=BL_COLORS.get(bl, "#555"), label=bl,
                   s=30, alpha=0.8)
    ax.axhline(0, ls=":", c="#444")
    rho, p = stats.spearmanr(m["cal_compat_mean"], m["log10_fold"])
    ax.set_xlabel("Calibrated pollen compatibility")
    ax.set_ylabel("Seed yield for size (log10 fold of expectation)")
    ax.set_title(f"D. Seed yield vs calibrated compatibility (ρ = {rho:.2f}, p = {p:.2f})",
                 loc="left", weight="bold")
    ax.legend(fontsize=8)
    for ax in axes.ravel():
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Drift calibration on the SRK-genotyped sites — best fit Ne = {a:.3g} · N^{b:.1f}",
                 weight="bold")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"step31_drift_calibration.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
