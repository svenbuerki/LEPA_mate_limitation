#!/usr/bin/env bash
# Run the LEPA mate-limitation pipeline end to end, in dependency order.
# See PIPELINE.md for what each step reads and writes.
#
# Usage (from the repository root):
#   ./run_pipeline.sh                 # all parts
#   ./run_pipeline.sh units           # Part 1 only (populations, demes, BLs)
#   ./run_pipeline.sh design          # model-population sampling plan only
#   ./run_pipeline.sh talkmaps        # presentation maps only (R)
#   ./run_pipeline.sh bulking         # seed-bulking table, figures and maps
#   PYTHON=/path/to/python ./run_pipeline.sh
#
# External inputs are set in Scripts/project_paths.py (override with the
# LEPA_DB, SRK_TABLES or LEPA_PROJECTS_ROOT environment variables).
set -euo pipefail

cd "$(dirname "$0")"
PY="${PYTHON:-/Users/sven/anaconda3/bin/python}"
PART="${1:-all}"
YEARS=(2026 2025)   # 2025 last: the single-year sampling tables it leaves
                    # behind are the ones Parts 3b and 4 read.
mkdir -p logs

run() {
  local step="$1"; shift
  echo "▶ ${step} $*"
  "$PY" "Scripts/${step}.py" "$@" > "logs/${step}${1:+_$2}.log" 2>&1 \
    || { echo "✖ ${step} failed — see logs/${step}${1:+_$2}.log"; exit 1; }
}

if [[ "$PART" == all || "$PART" == units ]]; then
  # Part 1 — analysis units (populations, demes, Bottleneck Lineages)
  run step01_populations
  run units_cartoon                 # schematic: slick spot → deme → population
  run step02_year_predictions
  run step03_BL_clustering
  run step04_BL_landscape_metrics
  run step05_BL_definition          # final populationIDs + BLs
  run step06_populations_classified
  run step07_event_to_deme          # event → population → deme table
  run step08_population_maps
  run step09_BL_overview_map
  run step10_demes_per_population_year
  run step11_year_scatter
  run step12_predictions_by_BL_year
fi

if [[ "$PART" == all || "$PART" == seedyield ]]; then
  # Part 2 — seed-yield test: mate limitation vs SI escape
  run step13_size_seed_yield
fi

if [[ "$PART" == all || "$PART" == sampling ]]; then
  # Part 3a — sampling design, one field season at a time
  for y in "${YEARS[@]}"; do
    run step14_seed_sampling_per_mother --year "$y"
    run step15_event_location_sampling --year "$y"
    run step16_connectivity
    run step17_fragmentation_sampling --year "$y"
  done
  # Part 3b — sampling diagnostics (2025 tables)
  run step18_mating_pool_structure
  run step19_fragmentation_index
  run step20_pollinator_radius
  run step21_sampling_adequacy      # both years; drops undersampled populations
fi

if [[ "$PART" == all || "$PART" == validation ]]; then
  # Part 4 — location-level predictions + seed-count power (2025)
  run step22_predictions_core --year 2025
  # Part 5 — validation against adult SRK genotypes (Part C)
  run step24_EO_scale_validation
  run step25_partC_clean_overlap
  run step26_partC_hypothesis_decomposition
  run step27_radius_sweep_H1
  run step31_drift_calibration      # fits drift to the genotyped sites (reads step25, step02, step13)
fi

# Model-population sampling plan (deme selection + germplasm pick list)
if [[ "$PART" == all || "$PART" == design ]]; then
  run step28_model_population_design
  echo "▶ step28_model_population_design_map (R)"
  Rscript Scripts/step28_model_population_design_map.R > logs/step28_model_population_design_map.log 2>&1 \
    || { echo "✖ step28_model_population_design_map failed — see logs/step28_model_population_design_map.log"; exit 1; }
  run step23_power_analysis          # needs step13, step21 and step28 tables
fi

# Consistency checks — units, BLs, IDs, targets, code, tables (fails the run)
if [[ "$PART" == all ]]; then
  run step29_pipeline_checks
fi

# Presentation maps on the topographic background (R; downloads elevation
# tiles once, cached in map_cache/). Outputs stay local.
if [[ "$PART" == all || "$PART" == talkmaps ]]; then
  echo "▶ step30_talk_maps (R)"
  Rscript Scripts/step30_talk_maps.R > logs/step30_talk_maps.log 2>&1 \
    || { echo "✖ step30_talk_maps failed — see logs/step30_talk_maps.log"; exit 1; }
fi

# Seed bulking — reproductive status of every population, talk figures and
# maps (figures_seed_bulking/, local); needs the seed-yield and sampling parts
if [[ "$PART" == all || "$PART" == bulking ]]; then
  run step32_seed_bulking
  echo "▶ step32_seed_bulking_maps (R)"
  Rscript Scripts/step32_seed_bulking_maps.R > logs/step32_seed_bulking_maps.log 2>&1 \
    || { echo "✖ step32_seed_bulking_maps failed — see logs/step32_seed_bulking_maps.log"; exit 1; }
fi

echo "✔ done (logs in logs/)"
