# LEPA_mate_limitation

**SRK mate-limitation and spatial-genetic framework for *Lepidium papilliferum* (slickspot peppergrass).**

This repository turns the per-individual S-allele genotypes produced by the Canu-amplicon pipeline into a **spatial-genetic framework** (populations, demes, Bottleneck Lineages), **predictions of SRK allele diversity and pollen compatibility** at every population and every year, and a **conservation-priority classifier** that bins each population as *mate-limited*, *reproduction intact*, or *SI-escape candidate* using the behavioural seed-yield signal (no SRK genotype required). The working hypothesis: **seed production reflects the reproductive status of a population** — under-production against a size-based expectation signals mate limitation, over-production signals SI escape.

The core deliverable is a priority list of populations where mate limitation is surfacing *right now* — the populations that should be prioritised for seed-genotyping, cross-pollination intervention, or founder augmentation. The companion deliverable is the biological interpretation: **why some small populations collapse while others buffer** under the same drift pressure.

---

## Where this fits in the broader LEPA SRK programme

Three repositories together carry the LEPA SRK conservation-genomics pipeline. Each is linked here by purpose; they can be read in order:

1. **[`svenbuerki/SRK_bioinformatics`](https://github.com/svenbuerki/SRK_bioinformatics)** — Phases 1 – 4 of the pipeline. Oxford Nanopore amplicon reads → per-individual S-allele genotypes → per-individual self-incompatibility status (SI / pSI / SC) → forward-time inheritance projections → experimentally executable cross plan. The *inputs* this repo consumes (species-wide prior `P1`, functional-allele definitions, SI-status calls, Canu-amplicon zygosity distribution) all come from the Phase 4 outputs of that repo.
2. **This repo — `svenbuerki/LEPA_mate_limitation`** — spatial-genetic framework + mate-limitation predictions + behavioural validation. Builds **populations**, **Bottleneck Lineages**, and **operational demes** on top of the Phase 4 genotype calls, generates per-population per-year predictions of SRK allele diversity and pollen compatibility, and classifies every population as mate-limited vs buffered via a species-wide size → seed-yield allometry.
3. **[`svenbuerki/Genetic-Rescue-DB`](https://github.com/svenbuerki/Genetic-Rescue-DB)** — the LEPA field database (SQL) and genetic-rescue experimental protocols. Supplies the per-event occurrence records, per-plant phenotyping (crown + height), per-mother germplasm records, and seed-yield measurements this repo consumes. The behavioural seed-yield test in `step13_size_seed_yield.py` (`Figures 12a / 12b` of the compact doc) was adapted from the `size_seed_model.py` family of scripts that originally lived in that repository.

```
                     ┌──────────────────────────────────┐
                     │  Genetic-Rescue-DB (field + lab) │
                     │  Occurrences · phenotyping ·     │
                     │  seed counts · germplasm records │
                     └───────────────┬──────────────────┘
                                     │
                                     ▼
┌──────────────────────────────────┐                 ┌───────────────────────────┐
│  LEPA_mate_limitation            │  ───────────▶   │  SRK_bioinformatics       │
│  Populations · Demes · BLs ·     │                 │  Phases 1 – 4 + Phase 5   │
│  predictions of SRK diversity    │                 │  cross design             │
│  + pollen compatibility ·        │  ◀───────────   │  Nanopore amplicons →     │
│  behavioural seed-yield test     │                 │  per-individual genotypes │
└──────────────────────────────────┘                 └───────────────────────────┘
```

**Reading of the two arrows:** the data flow between these two repos is **bidirectional** — neither side runs without the other.

- **`──▶` predictions flow out of this repo.** `LEPA_mate_limitation` produces, for every one of the 44 operational populations and every year, a predicted SRK allele diversity + a predicted pollen compatibility (sporophytic Class I / Class II dominance with empirical zygosity). These predictions feed `SRK_bioinformatics` Phase 5 — the hypothesis-driven cross design + validation step at the three populations that currently carry adult SRK genotypes (P6 = EO67, P40 = EO76, P42 = EO70).
- **`◀──` empirical inputs flow back in.** `SRK_bioinformatics` Phases 1 – 4 supply the ingredients the mate-limitation models are calibrated on: the species-wide prior `P1` (32 functional Fgs + their empirical frequencies), the functional-allele definitions, the Canu-amplicon zygosity distribution, and the per-individual SI / pSI / SC status calls. Change any of those inputs and this repo's predictions move.

The three repositories are intentionally **cross-linked by README reference only** — no git submodules, no DOI coupling. Each repository tracks its own publication timeline.

---

## What's in here

- **`Project_outline.md`** — the colleague-facing main document (**start here**): breeding units, predictions of SRK allele diversity and pollen compatibility, the seed-yield test and restoration priority list, the genotyped populations, model populations and discussion.
- **`Sampling_design.md`** — the agreed sampling plan for the D1 model populations (P36 + P34): why they are sampled, the two-step plan (mothers per deme, then spatially separated demes), the simulations that support it, deme maps, and the lab sheet of mothers to germinate. The earlier detailed version is `archive/Sampling_design_detailed.md` (not maintained).
- **`Seed_bulking.md`** — companion to the NCA Symposium talk: reproductive status of every population from seed set for plant size, the reproductively intact populations as candidate sources for seed bulking, and the talk figures (`figures_seed_bulking/`; the fine-scale per-BL maps stay local).
- **`archive/Sampling_and_prediction.md`** — the earlier long technical reference doc (formulas, output schemas, § A.5 – § C.6 derivations). **Archived 2026-10-17 and no longer maintained**; `Project_outline.md` takes precedence where they differ.
- **`SRK_mate_limitation_references.bib`** — citations for the metapopulation + Wright-neighbourhood + sporophytic-SI literature the pipeline leans on.
- **`PIPELINE.md`** — how to run the pipeline (`./run_pipeline.sh`), the run order, what every step reads and writes, where outputs are saved, and the files it needs from other repositories (LEPA database, SRK_bioinformatics).
- **`run_pipeline.sh`** — runs all 26 steps in dependency order.
- **`Scripts/`** — the 26 numbered pipeline steps (`step01_populations.py` through `step27_radius_sweep_H1.py`, numbered in run order) plus shared helpers (`srk_units.py` populations + demes, `project_paths.py` external inputs, `srk_bl_constants.py`, `srk_si_model.py`).
- **`tables/`** — all output TSV outputs: population crosswalk, trend-class classification, deme structure, radius-sensitivity sweep, predictions per (population, year), size-seed strata, conservation-priority list, and the Part C validation TSVs.
- **`figures/`** — all PNG + PDF outputs: 17 figures covering the populations overview, BL clustering + dendrogram, trend-class dynamic, deme-size distribution, SRK diversity + pollen-compatibility predictions, SI model schematic, Part C clean-overlap anchor, hypothesis decomposition, radius-sweep H1 test, size-seed calibration, and the per-population conservation-priority forest plot.

---

## Pipeline overview

**To run it, see [`PIPELINE.md`](PIPELINE.md).** The science runs in five stages:

1. **Spatial framework** — define the three nested units (populations · Bottleneck Lineages · operational demes) on the pooled 2025 + 2026 occurrence records.
2. **Across-year dynamic** — classify every population by its 2025 → 2026 trajectory (crash / stable / growth / ambiguous / single-year).
3. **Per-population predictions on all 44 populations** — predicted SRK allele diversity (per 50 m deme, unioned to the population) + predicted pollen compatibility (sporophytic Class I / Class II dominance with empirical zygosity, size-weighted mean across demes).
4. **SRK-based validation (Part C)** — at the six sites with SRK genotypes of plants grown from wild seed: EO67, EO70 and EO76 as whole sites, plus three populations of EO18 and EO27 placed through each plant's seed lot: observed vs predicted SRK diversity + observed vs predicted pollen compatibility + hypothesis decomposition + radius sweep.
5. **Behavioural prediction + validation on every population** — fit a species-wide plant-size → seed-yield allometry, compare each population's observed seed yield against the species-wide expectation, bin every population as *mate-limited* (under-producing → intervention priority), *reproduction intact* (at expectation), or *SI-escape candidate* (over-producing → seed-parent genotyping priority).

Scripts are numbered `step01` – `step27` in run order (renumbered 2026-10-17, step21 inserted 2026-10-19; old → new map in [`PIPELINE.md`](PIPELINE.md) § 7). Output TSVs and figures keep their historical `step28*` – `step30*` filename prefixes so that links from the docs and from `SRK_bioinformatics` stay stable.

---

## Status

- First commit: 2026-10-15. Repository seeded from the mate-limitation subset of the parent `SRK_bioinformatics` repo at commit `7e7ac29` (`6494fbc..7e7ac29` chain on `svenbuerki/SRK_bioinformatics`, pushed 2026-10-08).
- Behavioural seed-yield test (`step13_size_seed_yield.py`, Figures 12a + 12b of `Project_outline.md`) covers all 44 populations. Population-level restoration priority list (both years combined, calls that hold with and without the 2026 year effect): **6 mate-limited** · **1 switching** (over → under, P19) · **3 candidate mate-limited** · **7 SI-escape candidates** · **8 at expectation** · **9 year-sensitive** · **3 insufficient data**.
- Next test queued: within-BL stable small-vs-big pair D1 = P36 (EO25-B) × P34 (EO18-7 + EO18-8), 351 mothers with cleaned seed in the LEPA DB across 2025 + 2026.
- Scope: Snake River Plain populations only (44 populations). Jarbidge Foothills populations are out of scope here.

## License

GPL-3.0 — see [`LICENSE`](LICENSE).
