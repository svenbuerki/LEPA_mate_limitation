"""srk_units.py — the project's analysis units: population and deme.

Every analysis downstream of Step 01 – 04 works on these two units:

- **population** — 500 m-separated connected component on pooled
  2025 + 2026 events (Step 01, renumbered by BL in Step 04). Source of
  truth: `tables/step29a_population_crosswalk.tsv`.
- **deme** — 50 m connected component of one year's events inside one
  population (Step 01). Source of truth for deme sizes:
  `tables/step29a_demes_per_population_year.tsv`.

Step 01 does not store which deme each event belongs to, so
`load_event_units()` rebuilds that mapping from the crosswalk with the
same rule and the same event order as Step 01, then checks every deme
(n_events, N_fertile) against Step 01's deme table. A mismatch stops
the run: downstream results would otherwise sit on different units.

Run as a script to write the event → deme table:

    python Scripts/srk_units.py      # → tables/step29a_event_to_deme.tsv
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from step01_populations import DEME_RADIUS_M, connected_components

TABLES = Path("tables")
CROSSWALK_TSV = TABLES / "step29a_population_crosswalk.tsv"
DEMES_TSV = TABLES / "step29a_demes_per_population_year.tsv"
EVENT_TO_DEME_TSV = TABLES / "step29a_event_to_deme.tsv"

UNIT_COLS = ["event_year", "eventID", "populationID", "BL", "demeID",
             "deme_n_events", "deme_N_fertile",
             "locationID", "locationCode", "n_fertile"]


def load_event_units(year: int | None = None,
                     validate: bool = True) -> pd.DataFrame:
    """One row per (event_year, eventID) with its populationID, BL and
    demeID. `year` restricts to one field season."""
    cw = pd.read_csv(CROSSWALK_TSV, sep="\t")
    rows = []
    for (pop_id, yr), sub in cw.groupby(["populationID", "event_year"],
                                        sort=True):
        sub = sub.reset_index(drop=True)
        comps = connected_components(sub["lat"].to_numpy(),
                                     sub["lon"].to_numpy(), DEME_RADIUS_M)
        for deme_id, comp in enumerate(comps, start=1):
            part = sub.iloc[comp]
            for _, r in part.iterrows():
                rows.append(dict(
                    event_year=int(yr), eventID=int(r["eventID"]),
                    populationID=int(pop_id), BL=r["BL"], demeID=deme_id,
                    deme_n_events=len(comp),
                    deme_N_fertile=int(part["n_fertile"].sum()),
                    locationID=int(r["locationID"]),
                    locationCode=r["locationCode"],
                    n_fertile=int(r["n_fertile"]),
                ))
    units = pd.DataFrame(rows, columns=UNIT_COLS)
    if validate:
        units = _match_step01_deme_ids(units)
    if year is not None:
        units = units[units["event_year"] == int(year)].reset_index(drop=True)
    return units


def _match_step01_deme_ids(units: pd.DataFrame) -> pd.DataFrame:
    """Check the rebuilt partition against Step 01's deme table and adopt
    Step 01's demeIDs. Within each (population, year) the demes are
    matched on (n_events, N_fertile); demes tied on both are
    interchangeable for every size-based use."""
    demes = pd.read_csv(DEMES_TSV, sep="\t")
    key = ["deme_n_events", "deme_N_fertile"]
    mapping = []
    bad = []
    for (pop_id, yr), ref in demes.groupby(["populationID", "year"]):
        reb = (units[(units["populationID"] == pop_id)
                     & (units["event_year"] == yr)]
               .drop_duplicates("demeID"))
        ref_k = (ref.rename(columns={"n_events": "deme_n_events",
                                     "component_N_fertile": "deme_N_fertile"})
                 .sort_values(key + ["demeID"]))
        reb_k = reb.sort_values(key + ["demeID"])
        if (len(ref_k) != len(reb_k)
                or not (ref_k[key].to_numpy() == reb_k[key].to_numpy()).all()):
            bad.append((pop_id, yr))
            continue
        mapping += [dict(populationID=pop_id, event_year=yr,
                         demeID=int(a), demeID_step01=int(b))
                    for a, b in zip(reb_k["demeID"], ref_k["demeID"])]
    n_units = units.groupby(["populationID", "event_year"]).ngroups
    if bad or n_units != demes.groupby(["populationID", "year"]).ngroups:
        raise SystemExit(
            f"[srk_units] deme partition differs from {DEMES_TSV.name} "
            f"at (population, year) = {bad[:10]} — re-run Steps 01 – 04 "
            f"before downstream steps.")
    units = units.merge(pd.DataFrame(mapping),
                        on=["populationID", "event_year", "demeID"])
    return (units.drop(columns="demeID")
                 .rename(columns={"demeID_step01": "demeID"})[UNIT_COLS]
                 .sort_values(["event_year", "populationID", "demeID",
                               "eventID"])
                 .reset_index(drop=True))


def sampling_year(tables: Path = TABLES) -> int:
    """The single field season held by the current per-year sampling
    tables (Steps 14 – 17 are run one year at a time; 2025 last)."""
    ev = pd.read_csv(tables / "step28_events_spatial_neighborhood.tsv", sep="\t",
                     usecols=["event_year"])
    years = sorted(set(ev["event_year"].dropna().astype(int)))
    if len(years) != 1:
        raise SystemExit(f"[srk_units] sampling tables hold years {years}; "
                         "run Steps 14 – 17 with a single --year")
    return years[0]


def location_deme_inventory(year: int) -> pd.DataFrame:
    """Demes of one year grouped by DB location — the replacement for the
    retired `step29c_event_to_component_50m.tsv`. One row per event with
    `component_id_50m` = deme index within its location (0-based, ordered
    by populationID, demeID). Demes never span two locations (checked)."""
    u = load_event_units(year)
    if (u.groupby(["populationID", "demeID"])["locationID"].nunique() > 1).any():
        raise SystemExit("[srk_units] a deme spans two DB locations; "
                         "location-level inventory is undefined")
    deme_keys = (u.drop_duplicates(["locationID", "populationID", "demeID"])
                   .sort_values(["locationID", "populationID", "demeID"]))
    deme_keys["component_id_50m"] = deme_keys.groupby("locationID").cumcount()
    u = u.merge(deme_keys[["populationID", "demeID", "component_id_50m"]],
                on=["populationID", "demeID"])
    return (u.rename(columns={"n_fertile": "event_n_fertile",
                              "deme_N_fertile": "component_N_fertile",
                              "deme_n_events": "component_N_events"})
             [["locationID", "locationCode", "eventID", "component_id_50m",
               "event_n_fertile", "component_N_fertile", "component_N_events",
               "populationID", "demeID"]])


def eo_population_label(eo: str, year: int = 2025) -> str:
    """Current population(s) of an EO's locations in one season, e.g.
    'P6' for EO67 or 'P40 + P41' for EO76. Empty if the EO has no event
    that year."""
    u = load_event_units(year)
    codes = u["locationCode"].astype(str)
    sub = u[(codes == eo) | codes.str.startswith(eo + "-")
            | codes.str.match(rf"^{eo}[A-Z]")]
    return " + ".join(f"P{p}" for p in sorted(sub["populationID"].unique()))


def main() -> None:
    units = load_event_units()
    units.to_csv(EVENT_TO_DEME_TSV, sep="\t", index=False)
    n_demes = units.groupby(["event_year", "populationID", "demeID"]).ngroups
    print(f"[srk_units] Wrote {EVENT_TO_DEME_TSV} — {len(units)} events, "
          f"{units['populationID'].nunique()} populations, {n_demes} demes "
          f"(validated against {DEMES_TSV.name})")
    print(units.groupby("event_year")
               .agg(events=("eventID", "size"),
                    populations=("populationID", "nunique"),
                    demes=("demeID", lambda s: units.loc[s.index]
                           .groupby(["populationID", "demeID"]).ngroups))
               .to_string())


if __name__ == "__main__":
    main()
