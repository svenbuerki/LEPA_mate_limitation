# LEPA_mate_limitation

**Phase 5 of the SRK bioinformatics pipeline for *Lepidium papilliferum* (slickspot peppergrass).**

This repository turns the per-individual S-allele genotypes produced by the Canu-amplicon pipeline into a **spatial-genetic framework** (populations, demes, Bottleneck Lineages), **predictions of SRK allele diversity and pollen compatibility** at every population and every year, and a **conservation-priority classifier** that bins each population as *mate-limited* vs *buffered* using the behavioural seed-yield signal (no SRK genotype required).

The core deliverable is a priority list of populations where mate limitation is surfacing *right now* — the populations that should be prioritised for seed-genotyping, cross-pollination intervention, or founder augmentation. The companion deliverable is the biological interpretation: **why some small populations collapse while others buffer** under the same drift pressure.

---

## Where this fits in the broader LEPA SRK programme

Three repositories together carry the LEPA SRK conservation-genomics pipeline. Each is linked here by purpose; they can be read in order:

1. **[`svenbuerki/SRK_bioinformatics`](https://github.com/svenbuerki/SRK_bioinformatics)** — Phases 1 – 4 of the pipeline. Oxford Nanopore amplicon reads → per-individual S-allele genotypes → per-individual self-incompatibility status (SI / pSI / SC) → forward-time inheritance projections → experimentally executable cross plan. The *inputs* this repo consumes (species-wide prior `P1`, functional-allele definitions, SI-status calls, Canu-amplicon zygosity distribution) all come from the Phase 4 outputs of that repo.
2. **This repo — `svenbuerki/LEPA_mate_limitation` (Phase 5)** — spatial-genetic framework + mate-limitation predictions + behavioural validation. Builds **populations**, **Bottleneck Lineages**, and **operational demes** on top of the Phase 4 genotype calls, generates per-population per-year predictions of SRK allele diversity and pollen compatibility, and classifies every population as mate-limited vs buffered via a species-wide size → seed-yield allometry.
3. **[`svenbuerki/Genetic-Rescue-DB`](https://github.com/svenbuerki/Genetic-Rescue-DB)** — the LEPA field database (SQL) and genetic-rescue experimental protocols. Supplies the per-event occurrence records, per-plant phenotyping (crown + height), per-mother germplasm records, and seed-yield measurements that Phase 5 consumes. The behavioural seed-yield test in `step30i` (`Figures 12a / 12b` of the compact doc) was adapted from the `size_seed_model.py` family of scripts that originally lived in that repository.

```
                     ┌──────────────────────────────────┐
                     │  Genetic-Rescue-DB (field + lab) │
                     │  Occurrences · phenotyping ·     │
                     │  seed counts · germplasm records │
                     └───────────────┬──────────────────┘
                                     │
                                     ▼
┌───────────────────────────┐   ┌──────────────────────────────────┐
│  SRK_bioinformatics       │──▶│  LEPA_mate_limitation (Phase 5)  │
│  (Phases 1 – 4)           │   │  Populations · Demes · BLs ·     │
│  Nanopore amplicons →     │   │  SRK diversity + pollen-compat   │
│  per-individual genotypes │   │  predictions · behavioural seed- │
│  + SI / pSI / SC status   │   │  yield test · priority classifier│
└───────────────────────────┘   └──────────────────────────────────┘
```

The three repositories are intentionally **cross-linked by README reference only** — no git submodules, no DOI coupling. Each repository tracks its own publication timeline.

---

## What's in here

- **`Phase5_SRK_summary_for_colleagues.md`** — the compact colleague-facing summary (**start here**). Executive summary + Methodology + Result + interpretation (organised as 5 topic blocks: *Populations · Demes · SRK allele diversity · Pollen compatibility · Seed yield*) + Next test + conservation-priority classifier.
- **`Phase5_SRK_sampling_and_prediction.md`** — the long technical reference doc. Formulas, output schemas, design decisions, and the full § A.5 – § C.6 derivations.
- **`Phase5_references.bib`** — citations for the metapopulation + Wright-neighbourhood + sporophytic-SI literature the pipeline leans on.
- **`Scripts/`** — the 25 numbered pipeline steps (`step28_*` through `step30i_*`) plus 2 shared helpers (`srk_bl_constants.py`, `srk_si_model.py`). A clean 1 → 25 renumbering is planned as a follow-up commit — see "Pipeline overview" below for the current step naming.
- **`tables/`** — all Phase 5 TSV outputs: population crosswalk, trend-class classification, deme structure, radius-sensitivity sweep, predictions per (population, year), size-seed strata, conservation-priority list, and the Part C validation TSVs.
- **`figures/`** — all Phase 5 PNG + PDF outputs: 17 figures covering the populations overview, BL clustering + dendrogram, trend-class dynamic, deme-size distribution, SRK diversity + pollen-compatibility predictions, SI model schematic, Part C clean-overlap anchor, hypothesis decomposition, radius-sweep H1 test, size-seed calibration, and the per-population conservation-priority forest plot.

---

## Pipeline overview

Phase 5 runs in five stages (see § 30.1 – § 30.3 and § C.0.a – § C.0.c of the long doc for the detail):

1. **Spatial framework** — define the three nested units (populations · Bottleneck Lineages · operational demes) on the pooled 2025 + 2026 occurrence records.
2. **Across-year dynamic** — classify every population by its 2025 → 2026 trajectory (crash / stable / growth / ambiguous / single-year).
3. **Per-population predictions on all 44 populations** — predicted SRK allele diversity (per 50 m deme, unioned to the population) + predicted pollen compatibility (sporophytic Class I / Class II dominance with empirical zygosity, size-weighted mean across demes).
4. **SRK-based validation (Part C)** — at the three populations with Phase 4 adult SRK genotypes (P3 = EO67, P41 = EO70, P39 = EO76): observed vs predicted SRK diversity + observed vs predicted pollen compatibility + hypothesis decomposition + radius sweep.
5. **Behavioural prediction + validation on every population** — fit a species-wide plant-size → seed-yield allometry, compare each population's observed seed yield against the species-wide expectation, bin every population as *mate-limited* (seed shortfall → priority) vs *buffered* (seed yield at or above expectation).

The current `step*` names preserve the historical SRK-pipeline numbering (step28 through step30i); a clean 1 → 25 renumber is planned.

---

## Status

- First commit: 2026-10-15. Repository seeded from the Phase 5 subset of the parent `SRK_bioinformatics` repo at commit `7e7ac29` (`6494fbc..7e7ac29` chain on `svenbuerki/SRK_bioinformatics`, pushed 2026-10-08).
- Behavioural seed-yield test (Step 30i, Figures 12a + 12b of the compact doc) covers all 44 populations. Current tally: **8 mate-limited** · **4 candidate mate-limited** · **12 buffered** · **8 buffered-surplus** · **5 insufficient data**.
- Next test queued: within-BL stable small-vs-big pair D1 = P34 (EO25-B) × P32 (EO18-7 + EO18-8), 291 mothers available in the LEPA DB across 2025 + 2026.
- Scope: Snake River Plain populations only (39 Phase 5 locationCodes / 44 Phase 5 populations). Jarbidge Foothills populations are out of scope here.

## License

GPL-3.0 — see [`LICENSE`](LICENSE).
