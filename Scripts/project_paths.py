"""project_paths.py — single place for every input that lives OUTSIDE this repo.

All pipeline OUTPUTS go to `tables/` and `figures/` inside this repo (scripts
use paths relative to the repo root, so run them from there:
`python Scripts/stepNN_*.py`).

The three external inputs below are documented in PIPELINE.md (§ 5).
Each can be overridden with an environment variable, e.g.

    LEPA_DB=/path/to/LEPA_SQL.db python Scripts/step13_size_seed_yield.py
"""
from __future__ import annotations

import os
from pathlib import Path

# Parent folder holding the sibling projects on Sven's machine.
PROJECTS_ROOT = Path(os.environ.get(
    "LEPA_PROJECTS_ROOT", "/Users/sven/Documents/Current_projects"))

# 1. Production LEPA field + lab database (Genetic-Rescue-DB).
#    ESA-sensitive GPS data — local only, NEVER commit.
LEPA_DB = Path(os.environ.get(
    "LEPA_DB", PROJECTS_ROOT / "LEPA_fieldwork_protocol/SQL_DB/LEPA_SQL.db"))

# 2. SRK_bioinformatics git clone — Phase 4 outputs + EO-level BL tables.
SRK_TABLES = Path(os.environ.get(
    "SRK_TABLES",
    PROJECTS_ROOT / "SRK_bioinformatics/SRK_bioinformatics_GitHub/tables"))

PHASE4_INDIV_GENOTYPES_TSV = (
    SRK_TABLES / "Phase4/step23_individual_allele_genotypes_with_nulls.tsv")
EO_BL_SUMMARY_CSV = SRK_TABLES / "EO_BL_summary.csv"
EO_GROUP_BL_SUMMARY_CSV = SRK_TABLES / "EO_group_BL_summary.csv"
