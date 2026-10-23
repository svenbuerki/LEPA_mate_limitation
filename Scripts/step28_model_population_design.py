"""step28_model_population_design.py — deme-level sampling plan for the
model populations (default: D1 = P36 EO25-B small + P34 EO18-7/8 big).

For each population and field season:

1. **Fit demes.** A deme is dropped when it has a single fertile plant
   (no mate within 50 m, so no compatibility to measure) or when fewer
   usable mothers are held than its target (Step 21: mothers selected
   once every sampled seed lot is cleaned < deme target).
2. **Deme selection.** Each deme keeps its size-scaled target (Step 17:
   mothers for 90 % SRK allele detection, `srk_detection`). Demes are
   added greedily by expected new SRK alleles per mother, subject to a
   minimum separation between selected demes, until the expected share of
   the *whole* population's alleles detected (unfit demes included in the
   denominator) reaches the target.
   - Robustness: demes are chosen on one set of simulated deme pools and
     the stopping rule uses N_EVAL independent pool sets, so the reported
     coverage is not optimistic; the 5th percentile is reported too.
   - Stability: the selection is repeated on N_STABILITY pool sets and the
     share of runs that pick each deme is recorded.
   - Small populations: if all fit demes together need at most
     SMALL_POP_MOTHERS mothers, all are taken (spacing is moot).
3. **Pick list.** Within each selected deme, mothers are drawn at random
   (reproducible seed per population, year and deme) among those whose
   cleaned seed lot holds ≥ 25 seeds — one per slick spot, then the
   remaining places from the deme's other mothers; lots with fewer seeds,
   then uncleaned lots, are used only when eligible mothers run out. Seed
   count beyond 25 plays no part, so the sample is not biased towards the
   most fecund mothers. Each mother gets 25 seeds to germinate (≈ 15
   seedlings).

The model treats deme pools as independent draws from the species-wide
SRK frequencies; real demes have drifted further (Part C), so spatial
separation is a safeguard beyond what the model rewards.

Outputs (no coordinates): tables/step30_D_model_design_demes.tsv,
tables/step30_D_model_design_germplasm_by_deme.csv (lab list),
tables/step30_D_model_design_lab_sheet.tsv (lab sheet keyed on occurrenceID),
tables/step30_D_model_design_summary.tsv,
tables/step30_D_model_design_pick_list.tsv,
figures/step30_D_model_design_curves.{png,pdf}
"""
from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from project_paths import LEPA_DB
from srk_bl_constants import BL_COLORS
from srk_detection import (PLOIDY, PATERNAL_ALLELES_PER_SEEDLING,
                           SEEDLINGS_PER_MOTHER, species_prior)
from srk_units import CROSSWALK_TSV, load_event_units

TABLES = Path("tables")
FIGURES = Path("figures")

# ---- design settings (Sven, 2026-10-22) ----
MODEL_PAIRS = {"D1": {"small": 36, "big": 34},
               "D2": {"small": 21, "big": 19},
               "D3": {"small": 6, "big": 16}}
YEARS = (2025, 2026)
COVERAGE_TARGET = 0.90
MIN_SEPARATION_M = 100.0
MIN_PLANTS_PER_DEME = 2
SMALL_POP_MOTHERS = 10
SEEDS_TO_GERMINATE = 25
N_POOLS = 1000
N_EVAL = 20
N_STABILITY = 20
RNG_SEED = 2035


def deme_pools(sizes: np.ndarray, rng: np.random.Generator) -> list[np.ndarray]:
    """One simulated SRK pool per deme: (N_POOLS, K) allele frequencies
    from 4·N draws of the species prior."""
    f = species_prior()
    return [rng.multinomial(PLOIDY * int(n), f, size=N_POOLS) / (PLOIDY * int(n))
            for n in sizes]


def coverage(sel: list[int], pools: list[np.ndarray], A: np.ndarray) -> np.ndarray:
    """Per-pool share of the population's alleles (union over ALL demes)
    detected by the selected demes at their design effort A."""
    present = np.zeros(pools[0].shape, dtype=bool)
    for p in pools:
        present |= p > 0
    if not sel:
        return np.zeros(N_POOLS)
    miss = np.ones(pools[0].shape)
    for i in sel:
        miss *= (1.0 - pools[i]) ** A[i]
    return (1.0 - miss).sum(axis=1) / present.sum(axis=1)


def deme_distances(year: int, pop: int, deme_ids: list[int]) -> np.ndarray:
    """Shortest slick-spot-to-slick-spot distance (m) between demes."""
    cw = pd.read_csv(CROSSWALK_TSV, sep="\t", usecols=["event_year", "eventID", "lat", "lon"])
    u = load_event_units(year)
    ev = u[u["populationID"] == pop].merge(cw, on=["event_year", "eventID"])
    pts = {d: np.radians(ev.loc[ev["demeID"] == d, ["lat", "lon"]].to_numpy()) for d in deme_ids}
    n = len(deme_ids)
    D = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            p, q = pts[deme_ids[a]], pts[deme_ids[b]]
            dlat = q[None, :, 0] - p[:, None, 0]
            dlon = q[None, :, 1] - p[:, None, 1]
            h = (np.sin(dlat / 2) ** 2 + np.cos(p[:, None, 0]) * np.cos(q[None, :, 0])
                 * np.sin(dlon / 2) ** 2)
            D[a, b] = D[b, a] = (2 * 6_371_008.8 * np.arcsin(np.sqrt(h))).min()
    return D


def greedy(fit: list[int], pools_sel, pools_eval, A, M, D, min_sep):
    """Add demes by gain per mother (on pools_sel) until the mean coverage
    on the independent pools_eval reaches the target."""
    sel = []
    while True:
        ev = np.concatenate([coverage(sel, pe, A) for pe in pools_eval]) if sel else np.zeros(1)
        if ev.mean() >= COVERAGE_TARGET:
            break
        cand = [i for i in fit if i not in sel and all(D[i, j] >= min_sep for j in sel)]
        if not cand:
            break
        c0 = coverage(sel, pools_sel, A).mean()
        sel.append(max(cand, key=lambda i: (coverage(sel + [i], pools_sel, A).mean() - c0) / M[i]))
    return sel


def load_mothers() -> pd.DataFrame:
    con = sqlite3.connect(f"file:{LEPA_DB}?mode=ro", uri=True)
    df = pd.read_sql("""
        SELECT o.eventID, o.occurrenceID,
               MAX(g.germplasmID) AS germplasmID,
               COALESCE(MAX(g.germplasmQuantityEstimate), 0) AS seeds
        FROM Occurrences o JOIN Taxonomy t ON o.taxonID = t.taxonID
        LEFT JOIN Germplasm g ON g.occurrenceID = o.occurrenceID
        WHERE t.genus = 'Lepidium' AND t.specificEpithet = 'papilliferum'
          AND (o.provenance IS NULL OR o.provenance = 'in situ')
          AND o.eventID IS NOT NULL
        GROUP BY o.occurrenceID""", con)
    con.close()
    df["cleaned"] = df["germplasmID"].notna()
    germ = np.minimum(df["seeds"].astype(int), SEEDS_TO_GERMINATE)
    df["seedlings"] = np.where(df["cleaned"],
                               np.minimum(np.rint(germ * 0.60).astype(int), SEEDLINGS_PER_MOTHER),
                               SEEDLINGS_PER_MOTHER)
    return df[df["seedlings"] > 0]


def pick_mothers(m: pd.DataFrame, target: int,
                 rng: np.random.Generator) -> pd.DataFrame:
    """Random mothers within one deme, not ranked by seed yield.

    Tiers: (1) eligible — cleaned seed lot with ≥ SEEDS_TO_GERMINATE seeds
    (enough for the full 15 seedlings); (2) cleaned lot with fewer seeds
    (germinated in full); (3) seed not cleaned yet. One mother per slick
    spot is drawn at random from the best tier available in that slick
    spot; the remaining places are drawn at random from the deme's other
    mothers, best tier first. Seed count beyond the 25-seed threshold plays
    no part, so the sample is not biased towards the most fecund mothers.
    """
    m = m.copy()
    m["tier"] = np.select([m["cleaned"] & (m["seeds"] >= SEEDS_TO_GERMINATE), m["cleaned"]],
                          [1, 2], 3)
    m["draw"] = rng.random(len(m))
    m = m.sort_values(["tier", "draw"])
    floor = m.groupby("eventID", sort=False).head(1).assign(reason="one per slick spot (random)")
    rest = (m.drop(floor.index).head(max(target - len(floor), 0))
             .assign(reason="random among the deme's other mothers"))
    out = pd.concat([floor, rest])
    tier_note = {1: "", 2: "; lot < 25 seeds", 3: "; seed not cleaned yet"}
    out["reason"] = out["reason"] + out["tier"].map(tier_note)
    return out.drop(columns=["tier", "draw"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", default="D1", help="comma-separated model pairs (D1, D2, D3)")
    args = ap.parse_args()
    pairs = [p.strip() for p in args.pairs.split(",") if p.strip()]
    mothers = load_mothers()
    deme_rows, summary_rows, pick_rows, curves = [], [], [], []
    for pair in pairs:
        for role, pop in MODEL_PAIRS[pair].items():
            for year in YEARS:
                dt = pd.read_csv(TABLES / f"step30_A_sampling_adequacy_per_deme_{year}.tsv", sep="\t")
                d = dt[dt["populationID"] == pop].reset_index(drop=True)
                if d.empty:
                    continue
                d["fit"] = ((d["deme_N_fertile"] >= MIN_PLANTS_PER_DEME)
                            & (d["M_selected_capacity"] >= d["deme_M_target"]))
                d["unfit_reason"] = np.select(
                    [d["deme_N_fertile"] < MIN_PLANTS_PER_DEME,
                     d["M_selected_capacity"] < d["deme_M_target"]],
                    ["single plant", "too few mothers held"], "")
                M = d["deme_M_target"].to_numpy(int)
                A = M * (PLOIDY + PATERNAL_ALLELES_PER_SEEDLING * SEEDLINGS_PER_MOTHER)
                sizes = d["deme_N_fertile"].to_numpy(int)
                fit = [i for i in range(len(d)) if d.loc[i, "fit"]]
                D = deme_distances(year, pop, d["demeID"].astype(int).tolist())
                rng = np.random.default_rng([RNG_SEED, pop, year])
                pools_sel = deme_pools(sizes, rng)
                pools_eval = [deme_pools(sizes, rng) for _ in range(N_EVAL)]
                if M[fit].sum() <= SMALL_POP_MOTHERS:
                    sel, rule = list(fit), "all fit demes (small population)"
                else:
                    sel = greedy(fit, pools_sel, pools_eval, A, M, D, MIN_SEPARATION_M)
                    rule = f"greedy, ≥ {MIN_SEPARATION_M:.0f} m apart"
                cov_eval = np.concatenate([coverage(sel, pe, A) for pe in pools_eval])
                cov_all = np.concatenate([coverage(fit, pe, A) for pe in pools_eval])
                freq = np.zeros(len(d))
                if rule.startswith("greedy"):
                    for s in range(N_STABILITY):
                        r2 = np.random.default_rng([RNG_SEED, pop, year, s + 1])
                        ps = deme_pools(sizes, r2)
                        for i in greedy(fit, ps, [deme_pools(sizes, r2) for _ in range(3)],
                                        A, M, D, MIN_SEPARATION_M):
                            freq[i] += 1
                    freq /= N_STABILITY
                else:
                    freq[sel] = 1.0
                order = {i: k + 1 for k, i in enumerate(sel)}
                for i in range(len(d)):
                    r = d.loc[i]
                    deme_rows.append(dict(pair=pair, role=role, populationID=pop, year=year,
                                          demeID=int(r["demeID"]), fertile_plants=int(r["deme_N_fertile"]),
                                          mother_target=int(r["deme_M_target"]),
                                          mothers_held=int(r["n_occurrences"]),
                                          mothers_cleaned=int(r["n_cleaned"]),
                                          fit=bool(r["fit"]), unfit_reason=r["unfit_reason"],
                                          selected=i in order, selection_order=order.get(i),
                                          selected_in_share_of_runs=round(freq[i], 2)))
                cum, cm = [], 0
                for k, i in enumerate(sel):
                    cm += M[i]
                    cum.append((cm, float(np.concatenate([coverage(sel[:k + 1], pe, A)
                                                          for pe in pools_eval]).mean())))
                curves.append(dict(pair=pair, role=role, pop=pop, year=year, cum=cum,
                                   all_fit=(int(M[fit].sum()), float(cov_all.mean()))))
                summary_rows.append(dict(
                    pair=pair, role=role, populationID=pop, year=year, rule=rule,
                    demes=len(d), unfit_demes=int((~d["fit"]).sum()), fit_demes=len(fit),
                    mothers_all_fit=int(M[fit].sum()),
                    coverage_all_fit=round(float(cov_all.mean()), 3),
                    selected_demes=len(sel), mothers_selected=int(M[sel].sum()),
                    coverage_selected_mean=round(float(cov_eval.mean()), 3),
                    coverage_selected_p05=round(float(np.quantile(cov_eval, 0.05)), 3),
                    seeds_to_germinate=int(M[sel].sum()) * SEEDS_TO_GERMINATE,
                    seedlings_expected=int(M[sel].sum()) * SEEDLINGS_PER_MOTHER))
                # pick list
                u = load_event_units(year)
                mo = mothers.merge(u[["eventID", "populationID", "demeID"]], on="eventID")
                for i in sel:
                    dm = mo[(mo["populationID"] == pop) & (mo["demeID"] == d.loc[i, "demeID"])]
                    rng_pick = np.random.default_rng([RNG_SEED, pop, year,
                                                      int(d.loc[i, "demeID"]), 7])
                    for _, r in pick_mothers(dm, int(M[i]), rng_pick).iterrows():
                        pick_rows.append(dict(
                            pair=pair, role=role, populationID=pop, year=year,
                            demeID=int(d.loc[i, "demeID"]), eventID=int(r["eventID"]),
                            occurrenceID=int(r["occurrenceID"]),
                            germplasmID=(int(r["germplasmID"]) if r["cleaned"] else None),
                            status="germinate" if r["cleaned"] else "clean seed first",
                            seeds_available=int(r["seeds"]) if r["cleaned"] else None,
                            seeds_to_germinate=(min(int(r["seeds"]), SEEDS_TO_GERMINATE)
                                                if r["cleaned"] else SEEDS_TO_GERMINATE),
                            seedlings_expected=int(r["seedlings"]), reason=r["reason"]))
    demes = pd.DataFrame(deme_rows)
    summary = pd.DataFrame(summary_rows)
    picks = pd.DataFrame(pick_rows)
    picks["germplasmID"] = picks["germplasmID"].astype("Int64")
    picks["seeds_available"] = picks["seeds_available"].astype("Int64")
    demes.to_csv(TABLES / "step30_D_model_design_demes.tsv", sep="\t", index=False)
    summary.to_csv(TABLES / "step30_D_model_design_summary.tsv", sep="\t", index=False)
    picks.to_csv(TABLES / "step30_D_model_design_pick_list.tsv", sep="\t", index=False)
    # Lab list: the germplasmIDs to pull, per deme (CSV for the seed lab)
    lab = picks.sort_values(["pair", "role", "year", "demeID", "status", "germplasmID"])
    lab = lab.assign(population=lambda x: "P" + x["populationID"].astype(str))[
        ["pair", "role", "population", "year", "demeID", "germplasmID", "occurrenceID",
         "status", "seeds_available", "seeds_to_germinate", "seedlings_expected", "reason"]]
    lab.to_csv(TABLES / "step30_D_model_design_germplasm_by_deme.csv", index=False)
    # Lab sheet (TSV) keyed on occurrenceID — germplasmIDs are still missing
    # for seed lots not yet cleaned — with BL, population, EO and deme.
    u_all = load_event_units()
    ev_info = u_all.drop_duplicates("eventID")[["eventID", "BL", "locationCode"]]
    sheet = picks.merge(ev_info, on="eventID", how="left")
    sheet["EO"] = sheet["locationCode"].str.extract(r"^(EO\d+)")[0].fillna(sheet["locationCode"])
    sheet = sheet.assign(population="P" + sheet["populationID"].astype(str))[
        ["occurrenceID", "germplasmID", "BL", "population", "EO", "locationCode", "year",
         "demeID", "eventID", "pair", "role", "status", "seeds_available",
         "seeds_to_germinate", "seedlings_expected", "reason"]]
    sheet = sheet.rename(columns={"eventID": "eventID_slick_spot", "reason": "selected_because"})
    sheet = sheet.sort_values(["BL", "population", "year", "demeID", "eventID_slick_spot",
                               "occurrenceID"])
    sheet.to_csv(TABLES / "step30_D_model_design_lab_sheet.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print(f"[step28] pick list: {len(picks)} mothers "
          f"({(picks['status'] == 'germinate').sum()} ready, "
          f"{(picks['status'] != 'germinate').sum()} to clean first)")

    fig, axes = plt.subplots(1, len(curves), figsize=(3.6 * len(curves), 3.8), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, c in zip(axes, curves):
        uy = load_event_units(c["year"])
        col = BL_COLORS.get(uy.loc[uy["populationID"] == c["pop"], "BL"].iat[0], "#555")
        if c["cum"]:
            x, y = zip(*([(0, 0.0)] + c["cum"]))
            ax.plot(x, y, "-o", color=col, ms=4)
        ax.plot(*c["all_fit"], "s", color="#444", ms=6, mfc="white")
        ax.annotate("all fit demes", c["all_fit"], xytext=(4, -12),
                    textcoords="offset points", fontsize=7)
        ax.axhline(COVERAGE_TARGET, ls=":", c="#444", lw=0.8)
        ax.set_ylim(0, 1.02)
        ax.set_title(f"{c['pair']} {c['role']} — P{c['pop']}, {c['year']}", fontsize=9)
        ax.set_xlabel("Mothers sampled (cumulative)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Share of the population's SRK alleles\ndetected (independent simulations)")
    fig.suptitle("Model-population sampling: demes added one at a time (best gain per "
                 "mother, ≥ 100 m apart) until 90 % of the population's SRK alleles "
                 "are expected to be detected", fontsize=9)
    fig.tight_layout()
    FIGURES.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"step30_D_model_design_curves.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
