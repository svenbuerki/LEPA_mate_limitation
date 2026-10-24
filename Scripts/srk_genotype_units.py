"""srk_genotype_units.py — which population each SRK-genotyped plant belongs to.

Genotyped plants are grown ex situ from wild seed, or sampled as leaf
tissue from wild plants. The genotype table records only their EO; this
module places each plant in a project population (Step 01) through the
first link that resolves, in order:

1. **Seed lot** — the sample's germplasmID (master sequencing sheet) →
   `Germplasm` → its wild mother → her event → population.
2. **Location** — the sample's locationID (master sequencing sheet), when
   that location lies in a single population.
3. **Coordinates** — the tissue-bank position of the sampled plant,
   assigned to the population of the nearest event within 500 m (the
   population distance of Step 01).

Plants that resolve by none of these stay unplaced. Coordinates are read
and used in memory only; nothing with coordinates is written.
"""
from __future__ import annotations

import sqlite3

import numpy as np
import pandas as pd

from project_paths import (LEPA_DB, MASTER_SEQUENCING_XLSX,
                           PHASE4_INDIV_GENOTYPES_TSV, SAMPLING_METADATA_CSV,
                           TISSUE_BANK_XLSX)
from srk_units import CROSSWALK_TSV, load_event_units

POPULATION_RADIUS_M = 500.0
EARTH_R_M = 6_371_008.8


def _haversine(lat1, lon1, lat2, lon2):
    la1, lo1, la2, lo2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = (np.sin((la2 - la1) / 2) ** 2
         + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2)
    return 2 * EARTH_R_M * np.arcsin(np.sqrt(a))


def genotyped_plant_populations() -> pd.DataFrame:
    """One row per genotyped plant: Individual, EO, populationID (NaN if
    unplaced), `placed_by` (seed lot / location / coordinates /
    unplaced), and the seed lot (germplasmID) and occurrenceID it came
    from — seedlings sharing a seed lot are siblings of one wild mother."""
    geno = pd.read_csv(PHASE4_INDIV_GENOTYPES_TSV, sep="\t",
                       encoding="utf-8-sig", usecols=["Individual", "EO_normalised"])
    meta = pd.read_csv(SAMPLING_METADATA_CSV, encoding="utf-8-sig",
                       usecols=["SampleID", "OccurrenceID"])
    meta["OccurrenceID"] = pd.to_numeric(meta["OccurrenceID"], errors="coerce")
    x = (geno.merge(meta.drop_duplicates("SampleID"), left_on="Individual",
                    right_on="SampleID", how="left")
             .drop(columns="SampleID").rename(columns={"EO_normalised": "EO"}))

    master = pd.read_excel(MASTER_SEQUENCING_XLSX, sheet_name="Master Sheet",
                           usecols=["Occurrence", "Location", "Germplasm "])
    master["Occurrence"] = pd.to_numeric(master["Occurrence"], errors="coerce")
    master["germplasmID"] = pd.to_numeric(master["Germplasm "], errors="coerce")
    master["locationID"] = pd.to_numeric(master["Location"], errors="coerce")
    x = x.merge(master.dropna(subset=["Occurrence"])
                      .drop_duplicates("Occurrence")
                      [["Occurrence", "germplasmID", "locationID"]],
                left_on="OccurrenceID", right_on="Occurrence", how="left")

    units = load_event_units()
    ev_pop = units.drop_duplicates("eventID").set_index("eventID")["populationID"]
    loc_pop = units.groupby("locationID")["populationID"].agg(
        lambda s: s.iloc[0] if s.nunique() == 1 else np.nan)

    con = sqlite3.connect(f"file:{LEPA_DB}?mode=ro", uri=True)
    germ = pd.read_sql("""SELECT g.germplasmID,
                                 COALESCE(g.eventID, o.eventID) AS eventID
                          FROM Germplasm g
                          LEFT JOIN Occurrences o ON o.occurrenceID = g.occurrenceID""",
                       con)
    con.close()
    x["pop_seed_lot"] = x["germplasmID"].map(
        germ.set_index("germplasmID")["eventID"]).map(ev_pop)
    x["pop_location"] = x["locationID"].map(loc_pop)

    tissue = pd.read_excel(TISSUE_BANK_XLSX, sheet_name="Sheet1",
                           usecols=["OccurrenceID", "Latitude", "Longitude"])
    tissue = (tissue.dropna()
                    .assign(OccurrenceID=lambda d: d["OccurrenceID"].astype(int))
                    .drop_duplicates("OccurrenceID"))
    x = x.merge(tissue, on="OccurrenceID", how="left")
    cw = pd.read_csv(CROSSWALK_TSV, sep="\t", usecols=["lat", "lon", "populationID"])
    pop_xy = []
    for la, lo in zip(x["Latitude"], x["Longitude"]):
        if pd.isna(la) or pd.isna(lo):
            pop_xy.append(np.nan)
            continue
        d = _haversine(la, lo, cw["lat"].to_numpy(), cw["lon"].to_numpy())
        i = int(np.argmin(d))
        pop_xy.append(cw["populationID"].iat[i] if d[i] <= POPULATION_RADIUS_M else np.nan)
    x["pop_coordinates"] = pop_xy

    x["populationID"] = (x["pop_seed_lot"].fillna(x["pop_location"])
                                          .fillna(x["pop_coordinates"]))
    x["placed_by"] = np.select(
        [x["pop_seed_lot"].notna(), x["pop_location"].notna(),
         x["pop_coordinates"].notna()],
        ["seed lot", "location", "coordinates"], "unplaced")
    return (x[["Individual", "EO", "populationID", "placed_by", "germplasmID", "OccurrenceID"]]
              .drop_duplicates("Individual").reset_index(drop=True))
