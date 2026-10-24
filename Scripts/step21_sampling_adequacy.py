"""step21_sampling_adequacy.py — is each population sampled well enough to
test its predictions?

For every population and field season, take the mothers we actually hold
in each deme, select them the way Step 17 does (≥ 1 per event, then the
best-seeded up to the deme's `M_frag` target) and ask two questions:

1. **Coverage** — what fraction of the population's SRK allele classes
   (Fgs) would our seedlings detect? Each deme is simulated plant by
   plant (`srk_detection.py`: tetraploid with the empirical zygosity,
   alleles from the species prior); the deme's selected mothers show
   their own alleles plus 2 paternal alleles per genotyped seedling,
   fathers being compatible plants under sporophytic Class I / II SI.
   Population present / detected = union across demes — the model that
   also sets the design targets. 1 000 simulated populations → mean and
   95 % interval.
2. **Representation** — is **every deme** sampled (≥ 1 selected
   mother)? A deme with no mother is missing from every
   population-level estimate (the Fg union and the size-weighted
   P_compat mean), so the estimate is biased towards the sampled demes.
   The share of adults in unsampled demes is reported alongside.

Decision rule (agreed 2026-10-19): a population is **undersampled** and
dropped from the experiment for that season when mean coverage < 0.90
**or** any of its demes has no selected mother.

Mothers = LEPA in-situ occurrences (the plant we took seeds from). A
cleaned lot (germplasmID) gives min(15, round(0.6 × min(25, seeds)))
seedlings; an occurrence whose seeds are not cleaned yet is assumed to
give the design 15 (`capacity`). The `cleaned` columns count only
germplasm already in the DB, so the cleaning backlog stays visible.

Inputs:  LEPA DB (occurrences, germplasm) · tables/step29a_* (via
         srk_units) · tables/step29c_sampling_frag_aware_per_event_{Y}
         (deme targets, Step 17) · tables/step26i_L1_carrier_inventory
Outputs: tables/step30_A_sampling_adequacy_{Y}.tsv (one row per
         population) · tables/step30_A_sampling_adequacy_per_deme_{Y}.tsv
         · figures/step30_A_sampling_adequacy.{png,pdf}
"""
from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from project_paths import LEPA_DB
from srk_bl_constants import BL_COLORS
from srk_units import load_event_units
from step17_fragmentation_sampling import GERMINATION_RATE
from srk_detection import N_POOLS, deme_detection, popcount

TABLES = Path("tables")
FIGURES = Path("figures")
YEARS = (2025, 2026)
SEEDLINGS_PER_MOTHER = 15
SEEDS_TO_GERMINATE = int(np.ceil(SEEDLINGS_PER_MOTHER / GERMINATION_RATE))  # 25
COVERAGE_MIN = 0.90          # mean predicted Fg coverage
RNG_SEED = 2031


def load_mothers(db: Path = LEPA_DB) -> pd.DataFrame:
    """One row per LEPA in-situ occurrence: eventID, occurrenceID,
    cleaned (has a germplasmID) and seeds (germplasm estimate, 0 if not
    cleaned)."""
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    q = """
        SELECT o.eventID, o.occurrenceID,
               MAX(g.germplasmID IS NOT NULL)                AS cleaned,
               COALESCE(MAX(g.germplasmQuantityEstimate), 0) AS seeds
        FROM Occurrences o
        JOIN Taxonomy t ON o.taxonID = t.taxonID
        LEFT JOIN Germplasm g ON g.occurrenceID = o.occurrenceID
        WHERE t.genus = 'Lepidium' AND t.specificEpithet = 'papilliferum'
          AND (o.provenance IS NULL OR o.provenance = 'in situ')
          AND o.eventID IS NOT NULL
        GROUP BY o.occurrenceID"""
    df = pd.read_sql(q, con)
    con.close()
    return df


def seedlings_from_seeds(seeds: np.ndarray) -> np.ndarray:
    """Step 17's recipe: germinate min(25, seeds), keep ≤ 15 seedlings."""
    germ = np.minimum(seeds.astype(int), SEEDS_TO_GERMINATE)
    return np.minimum(np.rint(germ * GERMINATION_RATE).astype(int),
                      SEEDLINGS_PER_MOTHER)


def select_per_deme(mothers: pd.DataFrame, target: int) -> pd.DataFrame:
    """Step 17's selection inside one deme: the best-seeded mother of
    every event, then the best-seeded leftovers up to `target`."""
    m = mothers.sort_values(["seedlings", "occurrenceID"],
                            ascending=[False, True])
    floor = m.groupby("eventID", sort=False).head(1)
    rest = m.drop(floor.index).head(max(target - len(floor), 0))
    return pd.concat([floor, rest])


def deme_table(year: int, mothers: pd.DataFrame) -> pd.DataFrame:
    """One row per (population, deme) with adults, target and the
    selected mothers / seedlings under both mother definitions."""
    units = load_event_units(year)
    targets = pd.read_csv(
        TABLES / f"step29c_sampling_frag_aware_per_event_{year}.tsv",
        sep="\t")[["populationID", "demeID", "deme_N_fertile",
                   "deme_M_target"]].drop_duplicates(["populationID", "demeID"])
    demes = (units.drop_duplicates(["populationID", "demeID"])
                  [["populationID", "BL", "demeID", "deme_N_fertile"]]
                  .merge(targets, on=["populationID", "demeID"],
                         how="left", suffixes=("", "_step17")))
    bad = demes[demes["deme_N_fertile"] != demes["deme_N_fertile_step17"]]
    if len(bad):
        raise SystemExit(f"[step21] {year}: Step 17 deme table is out of date "
                         f"for {len(bad)} demes — rerun Step 17 --year {year}")
    demes["deme_M_target"] = demes["deme_M_target"].fillna(0).astype(int)

    mo = mothers.merge(units[["eventID", "populationID", "demeID"]],
                       on="eventID", how="inner")
    rows = []
    for _, d in demes.iterrows():
        sub = mo[(mo["populationID"] == d["populationID"])
                 & (mo["demeID"] == d["demeID"])].copy()
        row = d.drop("deme_N_fertile_step17").to_dict()
        row["n_occurrences"] = len(sub)
        row["n_cleaned"] = int(sub["cleaned"].sum())
        for scen in ("capacity", "cleaned"):
            s = sub.copy()
            if scen == "cleaned":
                s = s[s["cleaned"] == 1]
            s["seedlings"] = np.where(s["cleaned"] == 1,
                                      seedlings_from_seeds(s["seeds"].to_numpy()),
                                      SEEDLINGS_PER_MOTHER)
            s = s[s["seedlings"] > 0]
            sel = select_per_deme(s, int(d["deme_M_target"])) if len(s) else s
            row[f"M_selected_{scen}"] = len(sel)
            row[f"seedlings_{scen}"] = int(sel["seedlings"].sum()) if len(sel) else 0
        rows.append(row)
    out = pd.DataFrame(rows)
    out.insert(0, "year", year)
    return out


def simulate_population(demes: pd.DataFrame, scen: str,
                        rng: np.random.Generator) -> dict:
    """Fgs present (union of the demes' adults) and Fgs detected (union
    of what the selected mothers + seedlings recover) over simulated
    populations of `srk_detection` — the exact model the design targets
    come from, so a deme that meets its target scores what the design
    promised. Demes are simulated independently; each deme's mothers
    get its mean seedlings per mother. Intervals reflect drift and
    sampling uncertainty."""
    present = np.zeros(N_POOLS, np.uint64)
    detected = np.zeros(N_POOLS, np.uint64)
    for _, d in demes.iterrows():
        N = int(d["deme_N_fertile"])
        if N == 0:
            continue
        M = int(d[f"M_selected_{scen}"])
        s = int(round(d[f"seedlings_{scen}"] / M)) if M else 0
        p, det = deme_detection(N, M, N_POOLS, rng, seedlings=max(s, 1))
        present |= p
        detected |= det
    n_present = popcount(present)
    n_detected = popcount(detected)
    cov = np.divide(n_detected, n_present, out=np.zeros(N_POOLS),
                    where=n_present > 0)

    def q(a):
        return float(a.mean()), float(np.quantile(a, 0.025)), float(np.quantile(a, 0.975))
    out = {}
    for name, a in (("Fgs_present", n_present), ("Fgs_detected", n_detected),
                    ("coverage", cov)):
        out[f"{name}_mean_{scen}"], out[f"{name}_lo_{scen}"], out[f"{name}_hi_{scen}"] = q(a)
    return out


def classify(M: int, coverage: float, n_demes_unsampled: int) -> str:
    if M == 0:
        return "no_mothers"
    low_cov = coverage < COVERAGE_MIN
    low_rep = n_demes_unsampled > 0
    if low_cov and low_rep:
        return "undersampled_coverage+representation"
    if low_cov:
        return "undersampled_coverage"
    if low_rep:
        return "undersampled_representation"
    return "adequate"


def population_table(year: int, demes: pd.DataFrame) -> pd.DataFrame:
    units = load_event_units(year)
    codes = (units.groupby("populationID")["locationCode"]
                  .agg(lambda s: ", ".join(sorted(set(s)))))
    rows = []
    for pop_id, sub in demes.groupby("populationID"):
        N = int(sub["deme_N_fertile"].sum())
        row = dict(year=year, populationID=int(pop_id), BL=sub["BL"].iloc[0],
                   locationCodes=codes[pop_id], n_demes=len(sub), N_adults=N,
                   M_target=int(sub["deme_M_target"].sum()),
                   n_occurrences=int(sub["n_occurrences"].sum()),
                   n_cleaned=int(sub["n_cleaned"].sum()))
        rng = np.random.default_rng([RNG_SEED, year, int(pop_id)])
        for scen in ("capacity", "cleaned"):
            unsampled = sub[sub[f"M_selected_{scen}"] == 0]
            row[f"M_selected_{scen}"] = int(sub[f"M_selected_{scen}"].sum())
            row[f"seedlings_{scen}"] = int(sub[f"seedlings_{scen}"].sum())
            row[f"n_demes_unsampled_{scen}"] = len(unsampled)
            row[f"unsampled_adult_share_{scen}"] = (
                unsampled["deme_N_fertile"].sum() / N if N else 1.0)
            row.update(simulate_population(sub, scen, rng))
            row[f"status_{scen}"] = classify(
                row[f"M_selected_{scen}"], row[f"coverage_mean_{scen}"],
                row[f"n_demes_unsampled_{scen}"])
        row["include_in_experiment"] = row["status_capacity"] == "adequate"
        rows.append(row)
    return pd.DataFrame(rows)


def plot(pops: pd.DataFrame, out_png: Path, out_pdf: Path) -> None:
    """Rows = populations (both seasons), grouped BL1 → BL5. Panel A:
    predicted Fg coverage; B: adults in unsampled demes; C: mothers
    selected vs target. 2025 open, 2026 filled; red = undersampled."""
    labels = (pops.sort_values("year")
                  .groupby("populationID")
                  .agg(BL=("BL", "first"), codes=("locationCodes", "last")))
    labels = labels.sort_values(["BL", "populationID"])
    y_of = {p: i for i, p in enumerate(labels.index[::-1])}
    fig, axes = plt.subplots(1, 3, figsize=(15, 0.26 * len(labels) + 2.2),
                             sharey=True)
    style = {2025: dict(marker="o", mfc="white", offset=0.18),
             2026: dict(marker="o", mfc=None, offset=-0.18)}
    for yr, sub in pops.groupby("year"):
        st = style[yr]
        y = sub["populationID"].map(y_of).to_numpy() + st["offset"]
        ok = sub["include_in_experiment"].to_numpy()
        col = np.where(ok, sub["BL"].map(BL_COLORS).fillna("#555"), "#C00000")
        for ax, mean, lo, hi in (
                (axes[0], "coverage_mean_capacity", "coverage_lo_capacity",
                 "coverage_hi_capacity"),):
            ax.hlines(y, sub[lo], sub[hi], colors=col, lw=1, alpha=0.7)
        for i in range(len(sub)):
            mfc = "white" if yr == 2025 else col[i]
            r = sub.iloc[i]
            axes[0].plot(r["coverage_mean_capacity"], y[i], st["marker"],
                         mec=col[i], mfc=mfc, ms=5)
            axes[1].plot(r["unsampled_adult_share_capacity"], y[i], st["marker"],
                         mec=col[i], mfc=mfc, ms=5)
            axes[2].plot(r["M_target"], y[i], "|", color="#444", ms=7)
            axes[2].plot(r["M_selected_capacity"], y[i], st["marker"],
                         mec=col[i], mfc=mfc, ms=5)
            if r["M_selected_cleaned"] < r["M_selected_capacity"]:
                axes[2].plot(r["M_selected_cleaned"], y[i], "x",
                             color="#888", ms=4)
    axes[0].axvline(COVERAGE_MIN, ls=":", c="#444")
    axes[0].set_xlim(0, 1.02)
    axes[1].set_xlim(0, 1.0)
    axes[2].set_xscale("symlog", linthresh=10)
    axes[0].set_xlabel("Predicted share of the population's SRK alleles\ndetected by our seedlings")
    axes[1].set_xlabel("Share of fertile plants in demes\nwith no sampled mother")
    axes[2].set_xlabel("Mothers sampled (| = mothers needed;\nx = mothers with cleaned seed today)")
    for ax, t in zip(axes, ("A. SRK alleles detected", "B. Demes left unsampled",
                            "C. Mothers sampled vs needed")):
        ax.set_title(t, loc="left", weight="bold")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_yticks(range(len(labels)))
    axes[0].set_yticklabels([f"P{p} {labels.loc[p, 'codes']} ({labels.loc[p, 'BL']})"
                             for p in labels.index[::-1]], fontsize=7)
    for p in labels.index[::-1]:
        axes[0].get_yticklabels()[y_of[p]].set_color(
            BL_COLORS.get(labels.loc[p, "BL"], "#333"))
    fig.suptitle("Is each population sampled well enough? 2025 (open circles), "
                 "2026 (filled). Red = undersampled: our seedlings would detect "
                 "less than 90 % of its SRK alleles,\nor at least one of its "
                 "demes has no sampled mother. "
                 "BL = Bottleneck Lineage.", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    mothers = load_mothers()
    all_pops = []
    for yr in YEARS:
        demes = deme_table(yr, mothers)
        pops = population_table(yr, demes)
        demes.to_csv(TABLES / f"step30_A_sampling_adequacy_per_deme_{yr}.tsv",
                     sep="\t", index=False)
        pops.to_csv(TABLES / f"step30_A_sampling_adequacy_{yr}.tsv",
                    sep="\t", index=False)
        all_pops.append(pops)
        print(f"[step21] {yr}: {len(pops)} populations — "
              + ", ".join(f"{k} {v}" for k, v in
                          pops["status_capacity"].value_counts().items())
              + f" | cleaned today: "
              + ", ".join(f"{k} {v}" for k, v in
                          pops["status_cleaned"].value_counts().items()))
    FIGURES.mkdir(exist_ok=True)
    plot(pd.concat(all_pops), FIGURES / "step30_A_sampling_adequacy.png",
         FIGURES / "step30_A_sampling_adequacy.pdf")
    print("[step21] wrote tables/step30_A_sampling_adequacy_{2025,2026}.tsv, "
          "_per_deme_*, figures/step30_A_sampling_adequacy.{png,pdf}")


if __name__ == "__main__":
    main()
