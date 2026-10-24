#!/usr/bin/env python3
"""Step 29c — Fragmentation-aware sampling allocation.

Uses the **50 m primary pollen-flight radius** (see
`step29a_pollinator_radius_sensitivity.py` for the biological
justification) to group each location's events into connected
components (demes), then allocates mothers per deme (not per whole
location): the smallest M whose 15 seedlings per mother detect 90 % of
the deme's SRK Fgs on average, under the empirical species prior with
drift, tetraploid empirical zygosity and sporophytic Class I / II SI
(`srk_detection.py`; replaced the uniform-frequency coupon-collector
rule on 2026-10-19 and the independent-copy model on 2026-10-23). Within each component, mothers are
distributed proportional to N_fertile, subject to a maternal-
genotype floor of ≥ 1 mother per event.

Dependency order: run **Step 28 → Step 29 → Step 29b → Step 29c**.

The component-level total is summed to give the location-level
M_frag_aware, which is compared head-to-head with the existing
per-location recommendation in
step28_mothers_for_full_detection_by_location.tsv (§ B.4.1).

Why this is different
---------------------
The current allocation (Step 29 + B.4.1) uses the location's
*aggregate* N_fertile as a single pool. Two events sitting 300 m apart
at the same locationID are treated as one mating unit even though no
pollen crosses. The fragmentation-aware version scopes the coupon-
collector maths to each 50 m component, so:

  * A location that is one well-connected cluster → unchanged.
  * A location that is many isolated events → n_events mothers
    (private-allele floor was already forcing this).
  * A location that is a big cluster + a few isolated singletons
    (e.g. EO27-1, EO26-3) → the cluster only needs its coupon-
    collector M (~6), the isolated events each need 1, and the total
    can DROP substantially below the current recommendation while
    still delivering the same 90 %-see-every-allele guarantee.

Coupon-collector maths are scoped to each 50 m component rather than
the pooled location, so a location that is many isolated events is
allocated fewer mothers than the pooled Step 29 version would ask for,
and a location that is one well-connected cluster is unchanged.

Outputs
-------
    tables/step29c_sampling_frag_aware_per_event.tsv
    tables/step29c_sampling_comparison_per_location.tsv
    figures/step29c_sampling_comparison.png/pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from step14_seed_sampling_per_mother import (
    haversine_meters,
    mothers_for_full_detection,
    n_for_miss_probability,
    K_SPECIES_FG,
    K_pool,
    PLOIDY,
)
from srk_bl_constants import BL_COLORS, BL_ORDER, base_eo
from srk_units import load_event_units
from srk_detection import mothers_for_coverage

# Field → lab yield assumption used for Part C sampling design.
# Approximate LEPA germination rate under greenhouse conditions.
# Updating this constant propagates through n_seeds_to_germinate /
# n_seedlings_expected in step29c_partC_germplasmID_selection.tsv.
GERMINATION_RATE = 0.60

DEFAULT_TABLES = Path("tables")
DEFAULT_FIGURES = Path("figures")

EVENT_TSV = DEFAULT_TABLES / "step28_events_spatial_neighborhood.tsv"
R_PRIMARY = 50.0
TARGET_PROB = 0.90


def connected_components_50m(sub: pd.DataFrame) -> list[list[int]]:
    """BFS on the 50 m adjacency graph. Returns positional-index lists,
    one per component."""
    n = len(sub)
    if n == 0:
        return []
    lat = sub["lat"].values
    lon = sub["lon"].values
    d = haversine_meters(lat[:, None], lon[:, None],
                         lat[None, :], lon[None, :])
    adj = (d <= R_PRIMARY) & ~np.eye(len(d), dtype=bool)  # distance 0 links co-located events; only self is excluded
    seen = np.zeros(n, bool)
    comps: list[list[int]] = []
    for start in range(n):
        if seen[start]:
            continue
        stack = [start]
        comp: list[int] = []
        while stack:
            v = stack.pop()
            if seen[v]:
                continue
            seen[v] = True
            comp.append(v)
            for u in np.where(adj[v] & ~seen)[0]:
                stack.append(int(u))
        comps.append(comp)
    return comps


def allocate_component(comp_subset: pd.DataFrame,
                        target_prob: float = TARGET_PROB
                        ) -> pd.DataFrame:
    """Mothers the deme needs for 90 % mean Fg coverage under the
    empirical species prior with zygosity and SI
    (`srk_detection.mothers_for_coverage`), then distribute across events
    proportional to N_fertile with a floor of ≥ 1 per event and a cap
    at each event's N_fertile."""
    N_c = int(comp_subset["n_fertile"].sum())
    n_events_c = len(comp_subset)
    K_c = min(K_pool(N_c), K_SPECIES_FG) # reported only
    M_paternal = mothers_for_coverage(N_c) if N_c > 0 else 1
    M_c_final = max(M_paternal, n_events_c)

    events_out = comp_subset.copy().reset_index(drop=True)
    events_out["M_frag"] = 1 # maternal-genotype floor
    remaining = M_c_final - n_events_c
    if remaining > 0:
        weights = events_out["n_fertile"].astype(float).values
        weights = weights / weights.sum() if weights.sum() > 0 else weights
        raw = weights * remaining
        base = np.floor(raw).astype(int)
        residuals = raw - base
        extra = int(remaining - base.sum())
        if extra > 0:
            top = np.argsort(-residuals)[:extra]
            base[top] += 1
        events_out["M_frag"] = events_out["M_frag"].values + base

    # Cap at each event's N_fertile and redistribute any overflow
    n_cap = events_out["n_fertile"].astype(int).values
    while True:
        over_mask = events_out["M_frag"].astype(int).values > n_cap
        if not over_mask.any():
            break
        excess = int((events_out["M_frag"].values[over_mask]
                      - n_cap[over_mask]).sum())
        events_out.loc[over_mask, "M_frag"] = n_cap[over_mask]
        while excess > 0:
            headroom = (n_cap - events_out["M_frag"].astype(int).values)
            if (headroom > 0).sum() == 0:
                break
            i = int(np.argmax(headroom))
            events_out.at[i, "M_frag"] += 1
            excess -= 1
        if excess > 0:
            break

    events_out["component_K"] = K_c
    events_out["component_M_target"] = int(events_out["M_frag"].sum())
    events_out["component_N_events"] = n_events_c
    events_out["component_N_fertile"] = N_c
    events_out["M_paternal_only"] = M_paternal
    return events_out


def build_per_event(year: int) -> pd.DataFrame:
    """Allocate mothers per DEME (Step 01 50 m within-year component inside
    a population) — the drift unit and the pollen pool. Events come from
    Step 17's per-event table for `year`; their population and deme come
    from `srk_units`, so every allocation sits on the project units."""
    ev = pd.read_csv(EVENT_TSV, sep="\t", encoding="utf-8-sig")
    years = set(ev["event_year"].dropna().astype(int))
    if years != {year}:
        raise SystemExit(
            f"{EVENT_TSV} holds years {sorted(years)}, expected [{year}] — "
            f"run step17 and step11 with --year {year} first.")
    units = load_event_units(year)[["eventID", "populationID", "BL",
                                    "demeID", "locationCode"]]
    ev = ev.merge(units, on="eventID", how="left")
    orphans = ev[ev["populationID"].isna()]
    if len(orphans):
        print(f"[step29c] WARNING: {len(orphans)} {year} events are not in "
              f"any population (not in the Step 01 crosswalk) and are "
              f"skipped: eventIDs {sorted(orphans['eventID'].astype(int))}")
    ev = ev.dropna(subset=["populationID"]).copy()
    ev[["populationID", "demeID"]] = ev[["populationID", "demeID"]].astype(int)
    parts = []
    for (pop_id, deme_id), sub in ev.groupby(["populationID", "demeID"]):
        parts.append(allocate_component(sub.reset_index(drop=True)))
    return pd.concat(parts, ignore_index=True)


def build_per_population(event_df: pd.DataFrame,
                         partC: pd.DataFrame) -> pd.DataFrame:
    """Roll the deme allocation up to the population — the unit every
    result is reported on. `shortage` is the sum of DEME-level gaps: a
    mother from one deme cannot stand in for another deme's pollen pool."""
    deme = (event_df.drop_duplicates(["populationID", "demeID"])
            .groupby("populationID")
            .agg(n_demes=("demeID", "nunique"),
                 sum_M_paternal_only=("M_paternal_only", "sum"),
                 largest_deme_N_fertile=("component_N_fertile", "max"))
            .reset_index())
    pop = (event_df.groupby(["populationID", "BL"])
           .agg(locationCodes=("locationCode",
                               lambda s: ", ".join(sorted(set(map(str, s))))),
                n_events=("eventID", "count"),
                total_N_fertile=("n_fertile", "sum"),
                M_target=("M_frag", "sum"))
           .reset_index()
           .merge(deme, on="populationID"))
    pop["largest_deme_share"] = (pop["largest_deme_N_fertile"]
                                 / pop["total_N_fertile"]).round(3)
    avail = (partC.drop_duplicates(["populationID", "demeID"])
             .groupby("populationID")
             .agg(n_mothers_available_in_DB=("component_n_available_in_DB",
                                             "sum"),
                  shortage=("component_gap", "sum"))
             .reset_index())
    n_sel = partC.groupby("populationID").size().rename("n_mothers_selected")
    pop = (pop.merge(avail, on="populationID", how="left")
              .merge(n_sel, on="populationID", how="left"))
    for c in ("n_mothers_available_in_DB", "shortage", "n_mothers_selected"):
        pop[c] = pop[c].fillna(0).astype(int)
    # Demes with no germplasm at all never enter partC — their whole
    # target is a shortage.
    no_db = (event_df.drop_duplicates(["populationID", "demeID"])
             .merge(partC[["populationID", "demeID"]].drop_duplicates(),
                    on=["populationID", "demeID"], how="left", indicator=True))
    no_db = (no_db[no_db["_merge"] == "left_only"]
             .groupby("populationID")["component_M_target"].sum())
    pop["shortage"] = pop["shortage"] + pop["populationID"].map(no_db).fillna(0).astype(int)
    return pop.sort_values(["BL", "populationID"]).reset_index(drop=True)


def select_germplasm_for_partC(per_event_df: pd.DataFrame,
                                field_recipe_tsv: Path,
                                seeds_per_mother_target: int = 15
                                ) -> pd.DataFrame:
    """Pick specific germplasmIDs from the LEPA DB to satisfy step29c's
    fragmentation-aware allocation for Part C testing — **at the 50 m
    component scale**, not per event.

    Rationale. Events inside the same 50 m connected component share the
    same pollen pool (that is the biological definition of connectivity
    at the primary pollinator radius). A mother sampled at event A
    therefore samples the same pool as a mother at event B when both
    belong to the same component. So the coupon-collector coverage
    guarantee for the component travels freely across its events, and
    a shortage at one event can be absorbed by picking extra mothers
    at a connected event in the same component. Only the ≥ 1-mother-
    per-event *maternal-genotype floor* is a strict per-event
    requirement, and even that is waived at events that hold no
    germplasm in the DB at all.

    Algorithm (per 50 m component):
      1. Collect every germplasmID in the DB across every event of the
         component.
      2. Enforce the maternal-genotype floor: at each event that has at
         least one germplasmID, reserve the top-seeded mother.
      3. Fill the remaining `component_M_target − (# events with
         germplasm)` slots from the leftover mothers in the component,
         sorted by `seeds_available` (descending) then germplasmID
         (ascending) for a deterministic tie-break.
      4. If the component's total germplasm is smaller than
         `component_M_target`, take every available mother and record
         the shortfall in `component_gap`.

    Returns one row per SELECTED germplasmID with columns:
        germplasmID, occurrenceID, eventID, locationID, locationCode,
        seeds_available, n_seeds_to_genotype (= min(15, seeds_available)),
        component_id_within_loc, component_M_target,
        component_n_available_in_DB, component_gap,
        selection_reason ('floor' or 'top_seeds'),
        selection_priority_within_component,
        event_M_frag (for reference — the per-event allocation upstream).
    """
    if not field_recipe_tsv.exists():
        raise SystemExit(
            f"Missing {field_recipe_tsv} — run step29 first to build "
            "the per-germplasmID recipe.")
    recipe = pd.read_csv(field_recipe_tsv, sep="\t", encoding="utf-8-sig")

    # Bring per-event M_frag AND component identity + component target
    # from step29c's per-event output.
    event_meta = per_event_df[[
        "eventID", "populationID", "BL", "demeID", "M_frag",
        "component_M_target",
    ]].copy()
    # Mothers whose event is not in a population (see build_per_event
    # warning) cannot be assigned a deme and are left out.
    merged = recipe.merge(event_meta, on="eventID", how="inner")

    # Deme-level DB availability (all germplasmIDs at all events of the
    # deme — they share one pollen pool).
    comp_avail = (merged.groupby(["populationID", "demeID"])
                        .size().rename("component_n_available_in_DB")
                        .reset_index())
    merged = merged.merge(comp_avail, on=["populationID", "demeID"],
                          how="left")

    picked_rows = []
    # Iterate over demes.
    for (pop_id, deme_id), grp in merged.groupby(["populationID", "demeID"]):
        target = int(grp["component_M_target"].iloc[0])
        n_avail = int(len(grp))
        # Sort component-wide by seeds_available DESC, then germplasmID ASC.
        grp = grp.sort_values(
            ["seeds_available", "germplasmID"],
            ascending=[False, True]).reset_index(drop=True)

        # Step 1 — floor: one top-seeded mother per event with germplasm.
        floor_ids = (grp.groupby("eventID", sort=False)
                        .head(1)["germplasmID"].tolist())
        floor_set = set(floor_ids)
        floor_picked = grp[grp["germplasmID"].isin(floor_set)].copy()
        floor_picked["selection_reason"] = "floor"

        # Step 2 — fill remaining coverage slots.
        remaining_target = max(target - len(floor_picked), 0)
        leftover = grp[~grp["germplasmID"].isin(floor_set)]
        extra_picked = leftover.head(remaining_target).copy()
        extra_picked["selection_reason"] = "top_seeds"

        picked = pd.concat([floor_picked, extra_picked], ignore_index=True)
        picked["component_gap"] = max(target - n_avail, 0)
        picked_rows.append(picked)

    if picked_rows:
        selected = pd.concat(picked_rows, ignore_index=True)
    else:
        selected = merged.iloc[0:0].copy()
        selected["selection_reason"] = pd.Series(dtype="object")
        selected["component_gap"] = pd.Series(dtype=int)

    # Rank within component by seeds_available desc for the output.
    selected = selected.sort_values(
        ["populationID", "demeID",
         "selection_reason", "seeds_available", "germplasmID"],
        ascending=[True, True, True, False, True]).reset_index(drop=True)
    selected["selection_priority_within_deme"] = (
        selected.groupby(["populationID", "demeID"]).cumcount() + 1)

    # Part C genotyping is on SEEDLINGS, not seeds (user-confirmed
    # 2026-10-03). The tetraploid Rule 2 coupon-collector floor remains
    # 15 **seedlings** per mother; given a ~60 % germination rate we
    # need to germinate ~ceil(15 / 0.60) = 25 seeds per mother to hit
    # that seedling target in expectation. Columns:
    # seedlings_target — Rule 2 target (15 seedlings)
    # n_seeds_to_germinate — min(25, seeds_available)
    # n_seedlings_expected — round(n_seeds_to_germinate * 0.60)
    # n_seedlings_to_genotype — min(15, n_seedlings_expected)
    seeds_to_germinate_target = int(np.ceil(
        seeds_per_mother_target / GERMINATION_RATE)) # 25
    selected["seedlings_target"] = seeds_per_mother_target
    selected["n_seeds_to_germinate"] = np.minimum(
        selected["seeds_available"].astype(int),
        seeds_to_germinate_target,
    )
    selected["n_seedlings_expected"] = np.rint(
        selected["n_seeds_to_germinate"] * GERMINATION_RATE
    ).astype(int)
    selected["n_seedlings_to_genotype"] = np.minimum(
        selected["n_seedlings_expected"],
        seeds_per_mother_target,
    )
    # Backward-compat column (deprecated — kept for one release).
    selected["n_seeds_to_genotype"] = selected["n_seedlings_to_genotype"]

    # Project-wide EOID = base EO code (strips dash suffixes:
    # EO18-7 → EO18, EO27RT → EO27, EO27-3 → EO27). The TSV's primary
    # sort key is EOID; locationID disambiguates within-EO splits.
    selected["EOID"] = selected["locationCode"].astype(str).map(base_eo)

    keep = [
        "populationID", "BL", "demeID", "EOID", "locationCode", "locationID",
        "germplasmID", "occurrenceID",
        "eventID", "seeds_available",
        "seedlings_target", "n_seeds_to_germinate",
        "n_seedlings_expected", "n_seedlings_to_genotype",
        "n_seeds_to_genotype",
        "component_M_target", "component_n_available_in_DB",
        "component_gap", "selection_reason",
        "selection_priority_within_deme",
        "M_frag",
    ]
    keep = [c for c in keep if c in selected.columns]
    out = (selected[keep]
           .rename(columns={"M_frag": "event_M_frag"})
           .rename(columns={"component_M_target": "deme_M_target",
                            "component_n_available_in_DB":
                                "deme_n_available_in_DB",
                            "component_gap": "deme_gap"})
           .sort_values(["populationID", "demeID", "germplasmID"])
           .reset_index(drop=True))
    return out


def plot_comparison(loc_df: pd.DataFrame, out_png: Path, out_pdf: Path,
                    year: int):
    bls = [b for b in sorted(loc_df["BL"].dropna().unique()) if b in loc_df["BL"].values]
    if (loc_df["BL"] == "Unassigned").any():
        bls.append("Unassigned")
    palette = {**BL_COLORS, "Unassigned": "#8a8a8a"}
    heights = [max(int((loc_df["BL"] == b).sum()), 1) for b in bls]
    fig, axes = plt.subplots(
        len(bls), 1,
        figsize=(11.0, max(6.0, 0.36 * sum(heights) + 2.0)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )
    if len(bls) == 1:
        axes = [axes]
    for ax, bl in zip(axes, bls):
        sub = (loc_df[loc_df["BL"] == bl]
               .sort_values("M_target", ascending=True)
               .reset_index(drop=True))
        y = np.arange(len(sub))
        colour = palette[bl]
        h = 0.36
        # Solid bar — fragmentation-aware design target (M_frag).
        # Visually appears BELOW the light bar because barh with y - h/2
        # sits on the lower half of each row.
        ax.barh(y - h/2, sub["M_target"], h, color=colour, alpha=0.9,
                edgecolor=colour, linewidth=1.0)
        # Light bar — mothers already collected and stored in the LEPA DB.
        # Appears ABOVE the solid bar in each row.
        ax.barh(y + h/2, sub["n_mothers_available_in_DB"], h,
                color=colour, alpha=0.35, edgecolor=colour, linewidth=1.0)
        # Row label: populationID + deme structure.
        labels = [f"P{int(r['populationID'])} "
                  f"(events = {int(r['n_events'])}, "
                  f"demes = {int(r['n_demes'])}, "
                  f"adults = {int(r['total_N_fertile'])}, "
                  f"largest deme = {int(r['largest_deme_N_fertile'])})"
                  for _, r in sub.iterrows()]
        # Status annotation to the right of each row:
        # short by N → red, DB lacks mothers to hit the target
        # covered → grey, DB has enough (or more) already
        for i, r in sub.iterrows():
            short = int(r["shortage"])
            end_x = max(int(r["M_target"]),
                        int(r["n_mothers_available_in_DB"]))
            if short > 0:
                ax.text(end_x + 0.7, i, f"short by {short}",
                        fontsize=8, color="#c94b4b",
                        va="center", ha="left", fontweight="bold")
            else:
                ax.text(end_x + 0.7, i, "covered",
                        fontsize=8, color="#3c8f4c",
                        va="center", ha="left")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_ylim(-0.7, len(sub) - 0.3)
        ax.text(1.01, 0.5, bl, transform=ax.transAxes,
                fontsize=13, fontweight="bold", color=colour,
                va="center", ha="left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Figure-level legend with neutral grey handles so it doesn't look
    # BL-specific. Order matches the visual top-to-bottom layout in
    # each row: light 'available' bar on top, solid 'needed' bar below.
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="#555555", alpha=0.35, edgecolor="#333",
              label="Mothers already in the LEPA DB (light bar, on top)"),
        Patch(facecolor="#555555", alpha=0.9, edgecolor="#333",
              label="Mothers needed — 90 % SRK coverage per deme (solid bar, below)"),
    ]
    xmax = float(max(loc_df["M_target"].max(),
                     loc_df["n_mothers_available_in_DB"].max()))
    axes[-1].set_xlim(0, xmax + 18)
    axes[-1].set_xlabel(
        "Mothers per population — target (sum over its demes) vs already "
        "collected in the LEPA DB",
        fontsize=11)
    n_short_loc = int((loc_df["shortage"] > 0).sum())
    total_short = int(loc_df["shortage"].sum())
    total_target = int(loc_df["M_target"].sum())
    total_avail = int(loc_df["n_mothers_available_in_DB"].sum())
    fig.suptitle(
        f"{year} — mothers per population: do we have what the design asks for?\n"
        f"Target: {total_target}. In LEPA DB: {total_avail}. "
        f"Short: {n_short_loc} populations ({total_short} mothers, "
        f"counted per deme).",
        fontsize=11, y=0.998)
    fig.legend(handles=legend_handles, loc="upper center",
                bbox_to_anchor=(0.5, 0.94), ncol=2,
                fontsize=10, frameon=True)
    fig.tight_layout(rect=[0, 0, 0.94, 0.92])
    fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--year", type=int, required=True,
                    help="Field season (demes are within-year units). Run "
                         "step17 and step11 with the same --year first.")
    args = ap.parse_args()
    yr = args.year
    if not EVENT_TSV.exists():
        raise SystemExit(f"Missing input {EVENT_TSV}")

    ev = build_per_event(yr)
    ev_out = ev[[
        "populationID", "BL", "demeID", "locationID", "locationCode",
        "eventID", "n_fertile", "K_spatial_50m",
        "component_N_fertile", "component_N_events",
        "component_K", "component_M_target", "M_paternal_only",
        "M_frag"]].rename(columns={
            "component_N_fertile": "deme_N_fertile",
            "component_N_events": "deme_n_events",
            "component_K": "deme_K",
            "component_M_target": "deme_M_target"})
    out_ev = DEFAULT_TABLES / f"step29c_sampling_frag_aware_per_event_{yr}.tsv"
    ev_out.to_csv(out_ev, sep="\t", index=False)
    print(f"[step29c] Wrote {out_ev} ({len(ev_out)} events, "
          f"{ev_out.groupby(['populationID', 'demeID']).ngroups} demes, "
          f"{ev_out['populationID'].nunique()} populations)")

    # ---- Part C germplasmID selection (per deme) ----
    field_recipe_tsv = DEFAULT_TABLES / "step29_field_team_sampling_recipe.tsv"
    partC = select_germplasm_for_partC(ev, field_recipe_tsv)
    out_partC = DEFAULT_TABLES / f"step29c_partC_germplasmID_selection_{yr}.tsv"
    partC.to_csv(out_partC, sep="\t", index=False)

    pop = build_per_population(
        ev, partC.rename(columns={"deme_n_available_in_DB":
                                      "component_n_available_in_DB",
                                  "deme_gap": "component_gap"}))
    out_pop = DEFAULT_TABLES / f"step29c_sampling_per_population_{yr}.tsv"
    pop.to_csv(out_pop, sep="\t", index=False)

    n_target = int(pop["M_target"].sum())
    print(f"[step29c] {yr} sampling design (per deme, rolled up per population):")
    print(f" Target mothers                 : {n_target}")
    print(f" Mothers in LEPA DB (in demes)  : {int(pop['n_mothers_available_in_DB'].sum())}")
    print(f" Mothers selected               : {len(partC)}")
    print(f" Populations short of target    : {int((pop['shortage'] > 0).sum())} / {len(pop)}")
    print(f" Mothers short (sum over demes) : {int(pop['shortage'].sum())}")
    print(f" Seeds to germinate             : {int(partC['n_seeds_to_germinate'].sum())} "
          f"({GERMINATION_RATE:.0%} germination)")
    print(f" Seedlings to genotype          : {int(partC['n_seedlings_to_genotype'].sum())}")
    print(f"[step29c] Wrote {out_pop}, {out_partC}")

    DEFAULT_FIGURES.mkdir(parents=True, exist_ok=True)
    plot_comparison(
        pop,
        out_png=DEFAULT_FIGURES / f"step29c_sampling_per_population_{yr}.png",
        out_pdf=DEFAULT_FIGURES / f"step29c_sampling_per_population_{yr}.pdf",
        year=yr,
    )


if __name__ == "__main__":
    main()
