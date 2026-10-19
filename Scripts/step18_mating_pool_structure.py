"""Step 18 (historical tag: Step 29d) — mating-pool (deme) structure per
POPULATION, for the field season held by the sampling tables.

The framework defines a **50 m connected component** as one mating pool —
the set of adult plants that share a pollen pool at the 50 m
pollinator radius. A location may hold one mating pool (fully
connected) or many (fragmented), and each pool has its own
`component_N_fertile` adult count. This information is used by
step30's per-component simulations but has not had a dedicated
display figure.

This script fills that gap. It reads the project demes (Step 01 / 07,
via `srk_units`) for the sampling year and emits:

figures/step29d_mating_pool_structure.{png,pdf}
    Two-panel figure.
      Panel A — share of each population's adults in its largest
        deme (1.0 = one connected pool), one row per population.
      Panel B — per-population strip plot of deme sizes. One row per
        population, grouped vertically by BL in
        canonical BL_ORDER. Each dot = one 50 m mating pool, placed
        at its `component_N_fertile` on the log x-axis. Reference
        lines at N = 1 (single-plant SI floor) and N = 8 (coupon-
        collector 32-tetraploid-allele floor) — same convention as
        Figure 2b.

tables/step29d_mating_pool_summary.tsv
    One row per population with the deme structure (number of demes,
    largest deme, median deme size, etc.) summarised.

Dependency: Steps 01 – 07 (populations, demes) and the per-year sampling
tables (Steps 14 – 17; the year is read from them).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from srk_bl_constants import BL_COLORS, BL_ORDER
from srk_units import load_event_units, sampling_year

DEFAULT_TABLES = Path("tables")
DEFAULT_FIGURES = Path("figures")


def load_population_demes(year: int) -> pd.DataFrame:
    """Event-level deme table for one year, keyed for this script: the
    grouping unit (`locationID` column, kept for the plotting code) is the
    populationID, and `component_id_50m` is the demeID."""
    u = load_event_units(year)
    codes = (u.groupby("populationID")["locationCode"]
               .agg(lambda s: ", ".join(sorted(set(map(str, s))))))
    return pd.DataFrame({
        "locationID": u["populationID"],
        "locationCode": u["populationID"].map(lambda p: f"P{p}"),
        "locationCodes": u["populationID"].map(codes),
        "BL": u["BL"],
        "eventID": u["eventID"],
        "component_id_50m": u["demeID"],
        "component_N_fertile": u["deme_N_fertile"],
        "component_N_events": u["deme_n_events"],
    })


def build_summary(lookup: pd.DataFrame) -> pd.DataFrame:
    """Reduce the event-level lookup to one row per location, with
    mating-pool structure metrics."""
    # One row per (locationID, component_id_50m)
    pools = lookup.drop_duplicates(
        ["locationID", "component_id_50m"])[
        ["locationID", "locationCode", "locationCodes", "BL",
         "component_id_50m", "component_N_fertile", "component_N_events"]
    ].copy()
    pools["component_N_fertile"] = pools["component_N_fertile"].astype(int)

    rows = []
    for (loc_id, loc_code, codes, bl), sub in pools.groupby(
            ["locationID", "locationCode", "locationCodes", "BL"], sort=False):
        sizes = sub["component_N_fertile"].astype(int).values
        rows.append({
            "populationID": int(loc_id),
            "locationID": int(loc_id),
            "locationCode": loc_code,
            "locationCodes": codes,
            "BL": bl,
            "n_mating_pools_50m": int(len(sizes)),
            "total_adults": int(sizes.sum()),
            "largest_pool_N": int(sizes.max()),
            "smallest_pool_N": int(sizes.min()),
            "median_pool_N": float(np.median(sizes)),
            "largest_pool_share": float(sizes.max() / sizes.sum()),
            "n_pools_below_SI_floor_1": int((sizes <= 1).sum()),
            "n_pools_below_coupon_floor_8": int((sizes < 8).sum()),
        })
    summary = pd.DataFrame(rows)
    summary["BL"] = summary["BL"].fillna("Unassigned")
    summary["display_label"] = (summary["locationCode"] + " "
                                + summary["locationCodes"])

    return summary, pools


def _bl_order_for(summary: pd.DataFrame) -> list[str]:
    present = summary["BL"].unique().tolist()
    # Population-level BLs are numbered by rank (step05): BL1 → BL5.
    ordered = sorted((b for b in present if b != "Unassigned"),
                     key=lambda b: int(b[2:]))
    if "Unassigned" in present:
        ordered.append("Unassigned")
    return ordered


def plot_mating_pool_structure(summary: pd.DataFrame,
                                pools: pd.DataFrame,
                                out_png: Path, out_pdf: Path,
                                year: int) -> None:
    bls = _bl_order_for(summary)
    palette = {**BL_COLORS, "Unassigned": "#8a8a8a"}
    pools = pools.drop(columns="BL").merge(summary[["locationID", "BL"]],
                                           on="locationID", how="left")

    # Per-BL block heights proportional to the number of locations
    # in each BL.
    heights = [max(int((summary["BL"] == b).sum()), 1) for b in bls]
    fig = plt.figure(
        figsize=(15.0, max(6.5, 0.30 * sum(heights) + 1.8)),
    )
    # Outer gridspec: Panel A (connectivity share summary; carries
    # the row labels) + Panel B (per-location strip plot of pool
    # sizes; label-free, aligned with Panel A row for row).
    gs_outer = fig.add_gridspec(1, 2, width_ratios=[0.75, 1.4],
                                 wspace=0.06)
    gs_A = gs_outer[0, 0].subgridspec(len(bls), 1, hspace=0.14,
                                       height_ratios=heights)
    gs_B = gs_outer[0, 1].subgridspec(len(bls), 1, hspace=0.14,
                                       height_ratios=heights)
    axesA = [fig.add_subplot(gs_A[i, 0]) for i in range(len(bls))]
    axesB = [fig.add_subplot(gs_B[i, 0]) for i in range(len(bls))]

    # -----------------------------------------------------------------
    # Panel A — per-location connectivity share (= largest-pool share
    # of total adults). Single-number fragmentation diagnostic per
    # location: 1.0 = fully connected (whole census in one mating
    # pool), lower values = more fragmented.
    # -----------------------------------------------------------------
    # Pre-compute per-BL row ordering (sorted by largest_pool_N) so
    # both panels use the same order. Row labels live on Panel A.
    for ax_i, bl in zip(axesA, bls):
        sub = summary[summary["BL"] == bl].sort_values(
            "largest_pool_N", ascending=True).reset_index(drop=True)
        colour = palette[bl]
        y = np.arange(len(sub))
        share = sub["largest_pool_share"].astype(float).values
        ax_i.barh(y, share, color=colour, edgecolor="white",
                   height=0.72, zorder=2)
        ax_i.axvline(1.0, color="#444", ls=":", lw=0.9, alpha=0.6)
        ax_i.axvline(0.5, color="#b2182b", ls=":", lw=0.9, alpha=0.4)
        # Row labels on Panel A (leftmost column).
        y_labels = [
            f"{row['display_label']} "
            f"({row['n_mating_pools_50m']} deme"
            f"{'s' if row['n_mating_pools_50m'] != 1 else ''}, "
            f"{row['total_adults']} adults)"
            for _, row in sub.iterrows()
        ]
        ax_i.set_yticks(y)
        ax_i.set_yticklabels(y_labels, fontsize=8)
        ax_i.set_xlim(0.0, 1.05)
        ax_i.set_ylim(-0.6, len(sub) - 0.4)
        ax_i.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
        if ax_i is axesA[-1]:
            ax_i.set_xlabel(
                "Fraction in the largest deme "
                "(1.0 = fully connected)",
                fontsize=9,
            )
        else:
            ax_i.set_xticklabels([])
        ax_i.spines["top"].set_visible(False)
        ax_i.spines["right"].set_visible(False)
    axesA[0].set_title(
        "A Share of adults in the largest deme",
        fontsize=10, loc="left")

    # -----------------------------------------------------------------
    # Panel B — per-location strip plot: component_N_fertile per pool
    # -----------------------------------------------------------------
    for ax_i, bl in zip(axesB, bls):
        sub = summary[summary["BL"] == bl].sort_values(
            "largest_pool_N", ascending=True).reset_index(drop=True)
        colour = palette[bl]
        y_labels = []
        for y_pos, (_, row) in enumerate(sub.iterrows()):
            loc_pools = pools[pools["locationID"] == row["locationID"]]
            sizes = loc_pools["component_N_fertile"].astype(int).values
            # jitter dots vertically within the row for readability
            # Keep all dots on the row centre line — no vertical jitter,
            # dot size encodes pool size instead.
            ax_i.scatter(
                sizes, np.full(len(sizes), y_pos),
                s=np.clip(np.sqrt(sizes) * 4, 25, 180),
                color=colour, edgecolor="white", linewidth=0.5,
                alpha=0.85, zorder=3,
            )
            y_labels.append("") # row labels live on Panel A
        ax_i.axvline(1, color="#b2182b", ls=":", lw=1.0, alpha=0.7)
        ax_i.axvline(8, color="#333333", ls="--", lw=0.9, alpha=0.5)
        ax_i.set_xscale("log", base=2)
        ax_i.set_yticks(np.arange(len(sub)))
        ax_i.set_yticklabels(y_labels, fontsize=8)
        ax_i.set_ylim(-0.6, len(sub) - 0.4)
        ax_i.text(1.01, 0.5, bl, transform=ax_i.transAxes,
                   fontsize=12, fontweight="bold", color=colour,
                   va="center", ha="left")
        ax_i.spines["top"].set_visible(False)
        ax_i.spines["right"].set_visible(False)

    xticks = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
    axesB[-1].set_xticks(xticks)
    axesB[-1].set_xticklabels([str(x) for x in xticks], fontsize=9)
    axesB[-1].set_xlabel(
        "Adult plants per deme (log₂)",
        fontsize=10,
    )
    # Share x-limits
    all_sizes = pools["component_N_fertile"].astype(int)
    xmin = max(0.6, 0.9 * all_sizes.min())
    xmax = 1.15 * all_sizes.max()
    for ax_i in axesB:
        ax_i.set_xlim(xmin, xmax)
    # Only the last subpanel keeps x-tick labels; higher subpanels
    # share them visually.
    for ax_i in axesB[:-1]:
        ax_i.set_xticks(xticks)
        ax_i.set_xticklabels([])

    axesB[0].set_title(
        "B Deme sizes (one dot per 50 m component)",
        fontsize=10, loc="left")

    fig.suptitle(
        f"{year} — deme structure per population (50 m connected "
        "components within the year)",
        fontsize=12, y=0.995,
    )
    fig.subplots_adjust(left=0.22, right=0.96, top=0.93, bottom=0.08)
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    year = sampling_year()
    lookup = load_population_demes(year)
    summary, pools = build_summary(lookup)
    summary = summary.drop(columns="locationID")

    DEFAULT_TABLES.mkdir(parents=True, exist_ok=True)
    DEFAULT_FIGURES.mkdir(parents=True, exist_ok=True)

    tsv_path = DEFAULT_TABLES / "step29d_mating_pool_summary.tsv"
    summary.drop(columns=["locationCode"]).assign(year=year).to_csv(
        tsv_path, sep="\t", index=False)
    print(f"[step29d] Wrote {tsv_path} ({len(summary)} populations, "
          f"{len(pools)} demes, {year})")

    fig_png = DEFAULT_FIGURES / "step29d_mating_pool_structure.png"
    fig_pdf = DEFAULT_FIGURES / "step29d_mating_pool_structure.pdf"
    plot_mating_pool_structure(summary.assign(locationID=summary["populationID"]),
                               pools, fig_png, fig_pdf, year)
    print(f"[step29d] Wrote {fig_png} + .pdf")

    k = summary["n_mating_pools_50m"]
    s = pools["component_N_fertile"]
    print(f"[step29d] {year}: {len(summary)} populations, {len(pools)} demes "
          f"(mean {k.mean():.1f}, median {int(k.median())}, max {k.max()} per population); "
          f"deme sizes {s.min()} – {s.max()} (median {int(s.median())}); "
          f"single-plant demes {(s <= 1).sum()}; demes < 8 plants {(s < 8).sum()}")


if __name__ == "__main__":
    main()
