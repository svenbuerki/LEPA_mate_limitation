# Pipeline guide — LEPA mate limitation

How to run the pipeline, what each step reads and writes, where outputs
are saved, and which files come from other repositories and projects.
The read/write lists below were recorded by running the whole pipeline
under a file-access audit (2026-10-17), not compiled by hand.

---

## 1. Quick start

```bash
cd LEPA_mate_limitation           # always run from the repository root
./run_pipeline.sh                 # all 26 steps, in dependency order (~10 min)
./run_pipeline.sh units           # Part 1 only
./run_pipeline.sh seedyield       # Part 2 only
./run_pipeline.sh sampling        # Part 3 only
./run_pipeline.sh validation      # Parts 4 – 5 only
PYTHON=/path/to/python ./run_pipeline.sh
```

Each step writes a log to `logs/<step>.log` (not tracked by git). The run
stops at the first failing step.

**Python.** pandas, numpy, scipy, matplotlib, statsmodels, scikit-learn.
On Sven's machine the Anaconda interpreter has all of them and is the
script's default (`/Users/sven/anaconda3/bin/python`).

**Outputs.** Every table goes to `tables/` and every figure to
`figures/` (PNG + PDF), inside this repository.

**Coordinate tables stay local.** Nine tables carry precise coordinates
(`step28_events_spatial_neighborhood`, `step28_seed_sampling_per_mother`,
`step29a_population_crosswalk`, `step29a_slickspot_summary`,
`step30h_bl_definition`, `step30h_cluster_assignments`,
`step30h_landscape_per_BL`, `step30h_landscape_per_population`,
`step30h_populationID_crosswalk`). The pipeline writes them, but they are in
`.gitignore` because they give the locations of an ESA-listed plant. Never
force-add them. Figures, including maps, may be published.

---

## 2. Analysis units — everything depends on these

The pipeline defines its units once, in Part 1, and every later step
reads them.

| Unit | Definition | Built by | Source-of-truth table |
|---|---|---|---|
| **Event** | One occupied slick spot in one year (DB `Events`) | LEPA DB | — |
| **Slick spot** | Events matched across years within 10 m | step01 | `step29a_slickspot_summary.tsv` |
| **Population** | Connected component at **500 m** on pooled 2025 + 2026 events + historical locationID centroids | step01, numbered in step05 | `step29a_population_crosswalk.tsv` (event → population) |
| **Deme** | Connected component at **50 m** of one year's events inside one population (pollen pool, drift unit) | step01 | `step29a_demes_per_population_year.tsv` (sizes), `step29a_event_to_deme.tsv` (event → deme) |
| **Bottleneck Lineage (BL)** | Ward's D2 clustering of population centroids, *k* = 5, ordered by area then connectivity | step03 – step05 | `step30h_bl_definition.tsv`, `step30h_cluster_assignments.tsv` |

`populationID`s are assigned in step05: by BL (BL1 → BL5), then latitude
north → south within each BL. Re-running Part 1 after the data change can
renumber populations; `step30h_populationID_crosswalk.tsv` records
pre-numbering → final IDs for each run.

`Scripts/srk_units.py` (`load_event_units()`) is the shared accessor for
event → population → deme. It rebuilds the event → deme map from the
crosswalk with step01's rule and stops if it disagrees with step01's deme
table.

---

## 3. Run order and data flow

```
Part 1  Units        01 → 02 → 03 → 04 → 05 → 06 → 07 → (08 – 12 figures)
Part 2  Seed yield   13                                   needs 05, 06
Part 3  Sampling     for year in 2026, 2025: 14 → 15 → 16 → 17
                     then 18 → 19 → 20                    (2025 tables)
Part 4  Predictions  21 → 22                              (2025 tables)
Part 5  Validation   23 → 24 → 25 → 26                    needs 21
```

Ordering rules that are easy to get wrong:

- **step02 runs before step05.** step05 renumbers populations and
  rewrites the tables of steps 01 and 02 in place (it also backs them up
  to `tables/_bl_renumber_backup/`, not tracked).
- **step06 runs after step05** — it needs the final BL labels.
- **Part 3 runs per field season.** Demes are within-year units; pooling
  years merges plants that never share a pollen pool. step14 and step15
  take `--year`; step17 requires it. **2025 runs last**, because steps
  18 – 21 read the single-year `step28_*` / `step29_*` tables it leaves
  behind.
- **step21 runs after step17 for both years** (it reads both
  `step29c_sampling_frag_aware_per_event_{Y}` tables), and **step22
  runs after step21** (it takes the real mothers per deme from it).
- **step23 runs after step22** — it reads step22's
  `step30_A_prediction_prior_frequencies.tsv`.

Library imports are not run-order dependencies: steps 01 – 04, 15 – 17,
19, 20, 22 and 27 import helpers (haversine distance, DB loaders, coupon-
collector maths) from `step14_seed_sampling_per_mother.py`; steps 02,
21 and 27 import the drift model from `step22_predictions_core.py`;
steps 17, 20 and 21 size demes with `srk_detection.py` (mothers for
90 % mean Fg coverage under the empirical prior — the prediction model).

---

## 4. Steps — inputs and outputs

`t/` = `tables/`, `f/` = `figures/`. **EXT** = file outside this repo
(§ 5). Tables keep their historical `step28_` / `step29_` / `step30_`
filename prefixes so links from the docs and from `SRK_bioinformatics`
stay stable; § 7 maps them to the current step numbers.

### Part 1 — Analysis units (populations, demes, BLs)

| Step | Purpose | Reads | Writes |
|---|---|---|---|
| `step01_populations` | Slick spots, 500 m populations, 50 m within-year demes | EXT LEPA DB (`Events`, `Occurrences`, `Locations`, `Taxonomy`) · EXT EO BL tables | t/`step29a_population_crosswalk` · `step29a_population_summary` · `step29a_slickspot_summary` · `step29a_demes_per_population_year` · `step29a_seedcleaning_priority` |
| `step02_year_predictions` | Predicted SRK diversity + pollen compatibility per population per year (drift on demes) | t/`step29a_demes_per_population_year` · t/`srk_fg_class` · t/`srk_zygosity_empirical` · t/`step26i_L1_carrier_inventory` | t/`step30g_prediction_population_year` · `step30g_across_year_comparison` · `step30g_crash_candidates` · `step30g_stable_candidate_pair` |
| `step03_BL_clustering` | Ward's D2 clustering of population centroids, silhouette | t/`step29a_population_crosswalk` | t/`step30h_cluster_assignments` · `step30h_population_distances` · `step30h_silhouette_curve` · f/`step30h_dendrogram` · f/`step30h_silhouette_curve` |
| `step04_BL_landscape_metrics` | Area + connectivity per cluster (BL ordering keys) | t/`step29a_*` · t/`step30h_cluster_assignments` · t/`step30h_population_distances` | t/`step30h_landscape_per_population` · `step30h_landscape_per_BL` |
| `step05_BL_definition` | Final BL labels + populationIDs; rewrites all population tables | (the tables it rewrites) | rewrites 13 `step29a_*` / `step30g_*` / `step30h_*` tables · t/`step30h_bl_definition` · `step30h_populationID_crosswalk` |
| `step06_populations_classified` | Trend class per population (stable / crash / growth / single-year) | t/`step29a_population_crosswalk` · `step29a_population_summary` · `step30g_across_year_comparison` | t/`step30g_populations_classified` · f/`step30g_populations_classified` |
| `step07_event_to_deme` | Event → population → deme table (validated) | t/`step29a_population_crosswalk` · `step29a_demes_per_population_year` | t/`step29a_event_to_deme` |
| `step08_population_maps` | Population overview + zoom maps | t/`step29a_*` · `step30g_crash_candidates` · `step30g_stable_candidate_pair` | f/`step29a_populations_overview` · `step29a_crash_populations_zoom` · `step29a_candidate_populations_zoom` |
| `step09_BL_overview_map` | Populations coloured by BL | t/`step30h_bl_definition` · `step30h_cluster_assignments` | f/`step30h_overview_by_BL` |
| `step10_demes_per_population_year` | Deme sizes per population per year | t/`step29a_demes_per_population_year` | f/`step29a_demes_per_population_year` |
| `step11_year_scatter` | 2025 vs 2026 predictions per population | t/`step30g_across_year_comparison` | f/`step30g_across_year_scatter` |
| `step12_predictions_by_BL_year` | Predicted diversity + compatibility by BL and year | t/`step29a_population_crosswalk` · `step30g_prediction_population_year` | f/`step30h_pred_diversity` · `step30h_pred_pcompat` |

### Part 2 — Seed-yield test: mate limitation vs SI escape

| Step | Purpose | Reads | Writes |
|---|---|---|---|
| `step13_size_seed_yield` | Species-wide size → seed-yield curve (both years pooled); per population-year test; within-population year check — the test is repeated with that year effect as an offset on 2026, and statuses that differ are flagged `year_sensitive` (`robust_priority`); population-level restoration status; seed yield vs predicted compatibility / diversity | EXT LEPA DB (view `vOccurrenceTraits`) · t/`step29a_population_crosswalk` · `step30g_populations_classified` | t/`step30i_conservation_priority` (restoration priority list) · `step30i_size_seed_population_strata` · `step30i_size_seed_cv_rmse` · `step30i_year_effect_check` · f/`step30i_species_calibration` · `step30i_size_seed_population` · `step30i_seed_yield_vs_predictions[_tests]` · f/`step30i_seed_yield_vs_predictions` |

### Part 3 — Sampling design (per field season)

| Step | Purpose | Reads | Writes |
|---|---|---|---|
| `step14_seed_sampling_per_mother --year Y` | Seeds per mother (15 seedlings, 90 % allele detection); 50 m spatial neighbourhood per event | EXT LEPA DB (`Events`, `Occurrences`, `Germplasm`, `Taxonomy`) | t/`step28_seed_sampling_per_mother` · `step28_events_spatial_neighborhood` · `step28_coverage_curves_by_Nfertile` · `step28_aggregation_curves_by_Nfertile` · `step28_mothers_for_full_detection_by_bin` · `step28_mothers_for_full_detection_by_location` · f/`step28_coverage_curves` |
| `step15_event_location_sampling --year Y` | Per-mother field recipe (germplasm, seeds available) | EXT LEPA DB · t/`step28_seed_sampling_per_mother` · `step26i_L1_carrier_inventory` | t/`step29_field_team_sampling_recipe` · `step29_sampling_per_event` · `step29_sampling_per_location` · `step29_location_coverage_curves` |
| `step16_connectivity` | Connectivity at 10 / 25 / 50 m **per population** (sampling year) | t/`step28_events_spatial_neighborhood` · `step29a_event_to_deme` (via `srk_units`) | t/`step29_population_connectivity` |
| `step17_fragmentation_sampling --year Y` | **Mothers per deme**, rolled up per population; germplasmID selection | t/`step28_events_spatial_neighborhood` · `step29_field_team_sampling_recipe` · `step29a_population_crosswalk` · `step29a_demes_per_population_year` | t/`step29c_sampling_per_population_Y` · `step29c_sampling_frag_aware_per_event_Y` · `step29c_partC_germplasmID_selection_Y` · f/`step29c_sampling_per_population_Y` |
| `step18_mating_pool_structure` | Deme structure **per population** (sampling year) | t/`step29a_event_to_deme` (via `srk_units`) | t/`step29d_mating_pool_summary` · f/`step29d_mating_pool_structure` |
| `step19_fragmentation_index` | Fragmentation index per event / **population** | t/`step28_events_spatial_neighborhood` · `step29_population_connectivity` · `step29a_event_to_deme` | t/`step30_A_fragmentation_per_event` · `step30_A_fragmentation_per_population` · f/`step30_A_fragmentation_index` |
| `step20_pollinator_radius` | 10 – 200 m pollinator-radius sweep **per population**, both years (justifies 50 m) | EXT LEPA DB · t/`step29a_event_to_deme` · `step29c_sampling_per_population_2025\|_2026` · `srk_*` · `step26i_L1_carrier_inventory` | t/`step30_A_radius_sensitivity_summary[_2025\|_2026]` · `step30_A_radius_sensitivity_per_population_2025\|_2026` · f/`step30_A_radius_sensitivity` |
| `step21_sampling_adequacy` | Is each population sampled well enough? Real mothers per deme (Step 17 selection, empirical-prior targets) → predicted SRK Fg coverage + share of adults in unsampled demes; **undersampled** (coverage < 0.90 or any deme without a mother) populations are dropped, both years | EXT LEPA DB · t/`step29a_event_to_deme` · `step29c_sampling_frag_aware_per_event_2025\|_2026` · `step26i_L1_carrier_inventory` | t/`step30_A_sampling_adequacy_2025\|_2026` · `step30_A_sampling_adequacy_per_deme_2025\|_2026` · f/`step30_A_sampling_adequacy` |

### Part 4 — Location-level predictions + seed-count power (2025)

| Step | Purpose | Reads | Writes |
|---|---|---|---|
| `step22_predictions_core --year 2025` | Drift model per deme; predicted diversity + compatibility per **DB location** (Part C needs location level, see § 6); species prior | t/`step29_sampling_per_location` · `step29a_event_to_deme` (demes grouped by location, via `srk_units`) · `srk_*` · `step26i_L1_carrier_inventory` | t/`step30_A_prediction_location_diversity` · `step30_A_prediction_location_pcompat` · `step30_A_prediction_component_*` · `step30_A_prediction_prior_frequencies` · `step30_A_traffic_light_bands` · `step30_A_prediction_per_mother_fecundation` · f/`step30_A_prediction_fecundation` · `step30_A_diversity_vs_pcompat` · `step30_A_diversity_unbiased_vs_sampling` |
| `step23_seed_justification` | Precision + regression power at 15 seedlings per mother | t/`step30_A_prediction_prior_frequencies` | t/`step28d_pcompat_precision` · `step28d_matelim_power` · f/same names |

### Part 5 — Validation against adult SRK genotypes (Part C)

| Step | Purpose | Reads | Writes |
|---|---|---|---|
| `step24_EO_scale_validation` | Observed vs predicted compatibility at EO scale (split EOs) | EXT Phase 4 genotypes · `srk_*` · `step26i_L1_carrier_inventory` · `step30_A_prediction_prior_frequencies` | t/`step30_C_pcompat_validation_at_eo` · f/`step30_C_pcompat_observed_vs_predicted` |
| `step25_partC_clean_overlap` | Observed vs predicted at the three 1:1 populations (P3 EO67, P40 EO70, P38 EO76) | EXT Phase 4 genotypes · `step30_A_prediction_location_*` · `step30_A_traffic_light_bands` · prior | t/`step30_B_partC_clean_overlap_per_location` · `step30_B_partC_clean_overlap_fg_frequencies` · f/`step30_B_partC_clean_overlap` |
| `step26_partC_hypothesis_decomposition` | Which mechanism drives the compatibility gap | EXT Phase 4 genotypes · `step30_A_prediction_location_diversity` · bands · prior | t/`step30_B_partC_hypothesis_decomposition` · f/same name |
| `step27_radius_sweep_H1` | Can a tighter deme radius close the diversity gap? | EXT LEPA DB · t/`step30_B_partC_clean_overlap_per_location` · `step26i_L1_carrier_inventory` | t/`step30f_srk_diversity_radius_sweep` · f/same name |
| `step28_pipeline_checks` | Consistency checks; any failure stops the run: units rebuild = Step 01; BL names/order = Step 05; every table's populationID / deme / BL agrees with the units; deme targets = `srk_detection` rule; sampling tables = 2025; no legacy BL lookup or hard-coded BL order in code; every table has a writer or reader; no coordinate columns in tracked tables | all of tables/ · Scripts/ | console report (logs/step28_pipeline_checks.log) |

---

## 5. Dependencies on other repositories and projects

All external paths live in [`Scripts/project_paths.py`](Scripts/project_paths.py).
Each can be overridden with an environment variable, so the pipeline runs
on another machine without code edits:

```bash
LEPA_DB=/data/LEPA_SQL.db SRK_TABLES=/data/SRK_bioinformatics/tables ./run_pipeline.sh
LEPA_PROJECTS_ROOT=/data/projects ./run_pipeline.sh   # moves all defaults at once
```

### 5.1 Read live from outside this repository

| Input | Project / repository | Default path (under `/Users/sven/Documents/Current_projects/`) | Override | What is used | Steps |
|---|---|---|---|---|---|
| **LEPA production database** `LEPA_SQL.db` | [Genetic-Rescue-DB](https://github.com/svenbuerki/Genetic-Rescue-DB) (local `LEPA_fieldwork_protocol/`) | `LEPA_fieldwork_protocol/SQL_DB/LEPA_SQL.db` | `LEPA_DB` | `Events` (coordinates, `organismQuantityFertile`), `Occurrences`, `Locations`, `Taxonomy`, `Germplasm` (seed estimates), view `vOccurrenceTraits` (crown, height, seed yield) | 01, 13, 14, 15, 20, 26 |
| **Phase 4 per-individual SRK genotypes** `step23_individual_allele_genotypes_with_nulls.tsv` | [SRK_bioinformatics](https://github.com/svenbuerki/SRK_bioinformatics) | `SRK_bioinformatics/SRK_bioinformatics_GitHub/tables/Phase4/` | `SRK_TABLES` | Adult SRK alleles at EO67 / EO70 / EO76 and the split EOs | 23, 24, 25 |
| **EO-level BL tables** `EO_BL_summary.csv`, `EO_group_BL_summary.csv` | SRK_bioinformatics (originally `LEPA_EO_spatial_clustering`) | `SRK_bioinformatics/SRK_bioinformatics_GitHub/tables/` | `SRK_TABLES` | Legacy EO → BL lookup and BL colour order (`srk_bl_constants.py`); not used for the population-level BLs | most steps (via `srk_bl_constants`) |

> **The LEPA database holds ESA-sensitive GPS data. Never commit it or
> copy it into this repository** (`*.db` is in `.gitignore`).

Which events count: every event with coordinates and a fertile count
that either has a *L. papilliferum* in-situ occurrence or has no
occurrence records at all. The second case is an occupied slick spot whose
plants were recorded but not sampled (for example "plants too far gone
to collect"); its fertile plants are part of the deme and the census,
but it supplies no mothers or seed. step17 therefore cannot fill a deme
made only of such events (2026: P27 deme 2, P4 deme 2).

Data problems in the DB are reported, not guessed: events with a
non-numeric `organismQuantityFertile` are skipped by step14 with a
warning naming them; events without usable coordinates are left out of
every population by step01.

After a pipeline run that changes populations or demes, refresh the
database's `Populations`, `Demes`, `BottleneckLineages` tables and the
`Events.populationID` / `demeID` / `eventYear` columns with
`Genetic_Rescue_DB/Analyses/populations_and_demes/01_create_populations_demes_tables.py --apply`
(it reads `step29a_event_to_deme.tsv` and the population tables).

### 5.2 Copied into `tables/` from SRK_bioinformatics

These are static inputs. If they change upstream, copy the new version
into `tables/` and re-run from step02.

| File | Content | Upstream source | Steps |
|---|---|---|---|
| `srk_fg_class.tsv` | Class I / Class II assignment of the 32 functional SRK alleles (Fgs) | SRK_bioinformatics Phase 5 cross design | 02, 20, 21, 23 – 25 |
| `srk_zygosity_empirical.tsv` | Empirical tetraploid zygosity distribution | SRK_bioinformatics (Canu-amplicon genotypes) | 02, 20, 21, 23 – 25 |
| `step26i_L1_carrier_inventory.tsv` | Allele → functional-allele map (species-wide prior P1) | SRK_bioinformatics step26i | 02, 15, 20, 21, 23 – 26 |

### 5.3 What other projects read from here

`SRK_bioinformatics` Phase 5 (cross design) uses the per-population
predictions (`step30g_prediction_population_year.tsv`) and the
model-population choice. The Genetic-Rescue-DB `Populations`, `Demes` and
`BottleneckLineages` tables and the `Events.populationID` / `demeID`
columns mirror `step29a_event_to_deme.tsv`; refresh them after Part 1
changes.

---

## 6. Known limitations

- **Steps 21 and 23 – 26 (Part C) work by DB location / EO on purpose.**
  The adult SRK genotypes are recorded per EO only: the 885 genotyped
  occurrences in the LEPA DB have no eventID, locationID or date, so no
  adult can be placed in a population or deme. Step 22 therefore groups
  the current demes by location (demes never span two locations). Two of
  the three clean-overlap sites are one population (EO67 = P6, EO70 =
  P45); **EO76 is two populations (P43 + P44)** and is compared as one
  site. Linking genotyped occurrences to their events would allow a
  population-level Part C.
- All other steps work on populations and demes (`srk_units`); the
  retired `step29c_event_to_component_50m.tsv` is no longer read.
- **step14** labels locations from a previous step15 output when
  present (cosmetic; falls back to locationIDs on a first run).
- **step23** simulates the original 505-mother × 39-location design
  (hard-coded); power at the current 502-mother design is equivalent.

---

## 7. Script renumbering (2026-10-17)

Scripts were renumbered to match run order. Log lines and output
filenames keep the historical tags.

| Current | Previous | Historical tag (logs, output prefix) |
|---|---|---|
| step01_populations | step01 | step29a |
| step02_year_predictions | step14 | step30g |
| step03_BL_clustering | step03 | step30h (I) |
| step04_BL_landscape_metrics | step06 | step30h (II) |
| step05_BL_definition | step04 | step30h (III) |
| step06_populations_classified | step16 | step30g |
| step07_event_to_deme | — (new) | step29a |
| step08_population_maps | step02 | step29a |
| step09_BL_overview_map | step05 | step30h |
| step10_demes_per_population_year | step07 | step29a |
| step11_year_scatter | step15 | step30g |
| step12_predictions_by_BL_year | step20 | step30h |
| step13_size_seed_yield | step25 | step30i |
| step14_seed_sampling_per_mother | step17 | step28 |
| step15_event_location_sampling | step11 | step29 |
| step16_connectivity | step10 | step29b |
| step17_fragmentation_sampling | step12 | step29c |
| step18_mating_pool_structure | step08 | step29d |
| step19_fragmentation_index | step13 | step30_A |
| step20_pollinator_radius | step09 | step30_A |
| step22_predictions_core | step19 | step30_A / step30_B |
| step23_seed_justification | step18 | step28d |
| step24_EO_scale_validation | step21 | step30_C |
| step25_partC_clean_overlap | step22 | step30_B |
| step26_partC_hypothesis_decomposition | step23 | step30_B |
| step27_radius_sweep_H1 | step24 | step30f |
