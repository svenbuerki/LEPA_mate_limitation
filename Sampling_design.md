# Sampling design — testing the predictions with seed–mother genotyping

**Companion to [`Project_outline.md`](Project_outline.md).** The outline
defines the units (populations, Bottleneck Lineages, 50 m demes), the
predictions (SRK allele diversity and pollen compatibility) and the
seed-yield test. This document covers how we sample mothers and seeds to
test those predictions, and whether the material we already hold is
enough to do so. Run order and file details are in
[`PIPELINE.md`](PIPELINE.md).

---

## Summary — a two-step plan to recover 90 % of each population's SRK diversity

Testing the predictions needs **seed–mother genotyping**: each seedling
carries 2 paternal SRK alleles drawn from its mother's pollen pool, so the
seedlings of a mother sample the pollen pool of her deme. The goal is to
recover **90 % of the predicted SRK diversity of every population** while
sampling as few plants as possible and keeping a **balanced geographic
representation** of the population.

**Each mother: 25 seeds germinated → ≈ 15 seedlings genotyped** (§ 1).
Seeds germinate at ≈ 60 %; 15 seedlings detect every allele in a mother's
local pollen pool with 90 % probability and give the mate-limitation
regression 99.6 % power for a medium effect.

**Step 1 — within each deme: enough mothers to recover 90 % of the deme's
SRK alleles** (§ 2). The number of mothers is scaled to the deme's size
with the same drift model and species-wide allele frequencies as the
predictions: 2 mothers for a deme of 10 plants, 8 for 50, 13 for 100,
20 for 200 and 22 – 25 above 400, plus **at least one mother per slick
spot**. A deme is **unfit** — and set aside — if it has a single plant
(no mate within 50 m) or if we hold fewer usable mothers than it needs.

**Step 2 — across the population: spatially separated demes until 90 % of
the population's SRK alleles are recovered** (§ 4). Demes at least
**100 m apart** are added one at a time, the one adding the most expected
new alleles per mother first, until the selection is expected to detect
90 % of the whole population's SRK alleles. The spacing spreads the
sample across the population; the stopping rule keeps the number of
plants to the minimum. The coverage is judged on independent
simulations, so it is not optimistic.

**For the D1 model pair (P36 small + P34 big, BL3), the plan needs
93 mothers over 2025 + 2026 — 2 310 seeds to germinate, ≈ 1 390
seedlings** — instead of 191 mothers for every fit deme. 83 of the 93
have cleaned seed ready. Within each selected deme the mothers are
**drawn at random among those with enough seed (≥ 25)**, not ranked by
seed yield, so the sample represents the population's range of seed set
rather than its most fecund plants. The lab sheet, keyed on
occurrenceID, is
[`step30_D_model_design_lab_sheet.tsv`](tables/step30_D_model_design_lab_sheet.tsv).

**Can the material we hold test the predictions?** (§ 3) Across all
populations, 24 of 33 (2025) and 29 of 37 (2026) are sampled well enough
to recover 90 % of their SRK alleles with every deme represented; the
others are left out of that season's experiment.

---

## 1. Seedlings per mother

**Question.** How many seedlings must we genotype from each mother to
characterise her local pollen pool well enough to test the
mate-limitation prediction?

**Approach — per-mother allele detection.** Each seedling carries 2
paternal SRK alleles (tetraploid, sporophytic inheritance). Treating pollen
pickup as sampling with replacement from the local pollen pool, the chance
of missing one of K local alleles after 2·n seedling alleles is
`(1 − 1/K)^(2n)`. Setting it to 0.10 gives the seedlings needed for a 90 %
chance to see every allele. **This saturates at 15 seedlings per
mother** — beyond it, each extra seedling adds very little.

<a id="fig-s1"></a>
![Figure S1](figures/step28_coverage_curves.png)

**Figure S1.** SRK allele detection under tetraploid inheritance (4 SRK copies per plant; 2 paternal alleles per seedling), with the local pool capped at the species-wide 32 alleles. Shaded bands = 95 % simulation interval. **Panel A — one mother.** x = seedlings genotyped, one curve per slick-spot size class; the **red line at 15** marks the per-mother cap. One mother's 15 seedlings cannot saturate a 32-allele pool in large slick spots — the limit of a single mother, not undersampling. **Panel B — several mothers.** x = mothers sampled (15 seedlings each = 34 allele draws per mother: 4 maternal + 30 paternal). With 5 mothers, every size class reaches its local ceiling. Source: `step14_seed_sampling_per_mother.py`. Data: [`step28_coverage_curves_by_Nfertile.tsv`](tables/step28_coverage_curves_by_Nfertile.tsv), [`step28_aggregation_curves_by_Nfertile.tsv`](tables/step28_aggregation_curves_by_Nfertile.tsv).

**From seeds to seedlings — 60 % germination.** We genotype **seedlings**,
not seeds. LEPA seeds germinate at ≈ 60 % in the greenhouse, so reaching
15 genotyped seedlings per mother means germinating
`ceil(15 ÷ 0.60) = 25` seeds per mother. A seed lot with fewer than 25
seeds is germinated in full and delivers fewer seedlings.

**Is 15 seedlings enough statistical power?**

- **Per-mother pollen-compatibility precision.** Observed compatibility
  has binomial standard error `√(p·(1−p) / n)`, which falls from ≈ 0.28 at
  3 seedlings to ≈ 0.12 at 15 and then plateaus.
- **Mate-limitation regression power.** A full-pipeline simulation
  (505 mothers × 39 locations) gives, at 15 seedlings per mother,
  **99.6 % power for a medium effect** (100 seeds per unit compatibility)
  and 69 % for a small effect (50 seeds per unit, below the 80 % target).
  Errors-in-variables attenuation of the slope is ≈ 0.42 — real, but it
  does not prevent detection at realistic effect sizes.

<a id="fig-s2"></a>
![Figure S2](figures/step28d_matelim_power.png)

**Figure S2.** Statistical support for 15 seedlings per mother. **Left — per-mother pollen-compatibility precision:** standard error vs seedlings genotyped for four true compatibility values (0.3 / 0.5 / 0.7 / 0.9); precision gains plateau around 15. **Right — mate-limitation regression power** from a full-pipeline simulation, four effect sizes (50 / 100 / 200 / 400 seeds per unit compatibility); at 15 seedlings: 99.6 % power for the medium effect, 69 % for the small one. Source: `step14_seed_sampling_per_mother.py`. Companion precision figure: [`step28d_pcompat_precision.png`](figures/step28d_pcompat_precision.png).

---

## 2. Step 1 — mothers per deme

**Question.** How many mothers must we sample in a deme for their
seedlings to reveal the SRK alleles the deme holds?

**Approach.** The same model as the predictions (Project_outline § 30.1).
A deme of *N* fertile plants holds 4·*N* SRK alleles drawn from the
species-wide frequencies — drift, so rare alleles are often absent. *M*
mothers with 15 seedlings each sample 34·*M* alleles from that pool. The
target is the smallest *M* whose seedlings detect, on average over 2 000
simulated pools, **90 % of the SRK alleles present in the deme**. On top
of that, **at least one mother per slick spot**: an isolated slick spot's
private alleles cannot be recovered from any other slick spot.

An earlier rule treated all 32 alleles as equally common; it asked for
about 6 mothers in every deme of 10 or more plants — too many for small
demes and too few for large ones, where rare alleles need more draws:

| Deme size (fertile plants) | 10 | 20 | 50 | 100 | 200 | 420 | 1 000 |
|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| Mothers (species-wide frequencies, current) | 2 | 4 | 8 | 13 | 20 | 25 | 23 |
| Mothers (all alleles equally common, earlier) | 6 | 6 | 6 | 6 | 6 | 6 | 6 |
| Alleles detected with the earlier 6 mothers | 100 % | 97 % | 86 % | 77 % | 70 % | 66 % | 66 % |

**Totals across all populations.**

| Season | Populations / demes | Mothers needed | Mothers sampled (in DB) | Seed cleaned | Mothers usable | Short (counted per deme) | Seeds to germinate → seedlings |
|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 2025 | 33 / 99 | **695** | 765 | 765 | 549 | 146 | 13 625 → 8 176 |
| 2026 | 37 / 189 | **1 110** | 1 581 | 705 | 1 029 | 81 | ~25 700 → 15 339 |

A mother in one deme cannot stand in for another deme's pollen pool, so
shortfalls are counted deme by deme. 2026 counts every sampled mother; 876
of them still await seed cleaning, and their lots are assumed to deliver
15 seedlings each. Two 2026 demes (P25 deme 2, P4 deme 2) consist only of
slick spots recorded too late to collect seed: they are part of the model
but have no mothers to sample.

<a id="fig-s3"></a>
![Figure S3](figures/step29c_sampling_per_population_2025.png)

**Figure S3.** Do we hold the mothers the design asks for? 2025, one bar pair per population, grouped by Bottleneck Lineage (BL1 → BL5). **Solid bar** = mothers needed (sum over the population's demes); **light bar** = mothers with seed in the LEPA database. "Short by *n*" counts the shortfall deme by deme, so a population can hold more mothers than it needs in total and still be short in one deme (P34, P23). The 2026 version is [`step29c_sampling_per_population_2026.png`](figures/step29c_sampling_per_population_2026.png). Source: `step17_fragmentation_sampling.py`. Data: [`step29c_sampling_per_population_2025.tsv`](tables/step29c_sampling_per_population_2025.tsv).

---

## 3. Can our sampling test the predictions?

**Question.** With the mothers we actually hold, will the seedlings
represent each population — or is it undersampled and must be left out of
that season's experiment?

**Approach.** For each population and season, the mothers in each deme are
selected as the design does (at least one per slick spot, then the
best-seeded up to the deme's target), each with the seedlings her seed lot
can deliver (at most 15). The same drift model then predicts:

1. **What the population holds** — each deme's SRK alleles drawn from the
   species-wide frequencies (4 per plant); the population holds an allele
   if any of its demes does.
2. **What our sampling recovers** — in each deme, the mothers' 4 alleles
   plus 2 paternal alleles per seedling drawn from the deme's own allele
   frequencies; the population's recovered set is the union over demes.
3. **Coverage** = recovered ÷ held, averaged over simulated pools, with a
   95 % interval.

A deme with no mother adds its alleles to what the population holds but
nothing to what we recover — the sampling blind spot. **A population is
undersampled if its seedlings would recover less than 90 % of its SRK
alleles, or if any of its demes has no sampled mother**: the population is
represented only when every deme is sampled. Earlier versions split each
location's mothers across demes by size and counted every seed as
genotyped; both hid the blind spot.

**Result.**

| Season | Sampled well enough | Undersampled (left out) |
|---|:-:|---|
| 2025 | **24 of 33** | P2 (EO30-1), P8 (EO27), P10 (EO27-3), P13 (EO27-1) (BL1) · P26 (EO29, BL2) · P27 (EO48), P28, P29 (EO32) (BL3) · P41 (EO76, BL4) |
| 2026 | **29 of 37** | P2 (EO30-1), P4 (EO30-3), P10 (EO27-3), P13, P15 (EO27-1) (BL1) · P25 (EO61), P26 (EO29) (BL2) · P40 (EO76, BL4) |

- **Every deme has a sampled mother** in 2025 (99 of 99) and in 2026 except
  the two recorded-only demes (P4 deme 2, 62 plants; P25 deme 2, 20
  plants), where the crew arrived after seed had dropped.
- The other undersampled populations hold mothers in every deme but too
  few in at least one — mostly large single-deme or few-deme populations
  where the deme target exceeds the mothers collected (P26 EO29: 25 needed,
  22 held). Two mate-limited populations are among them — **P10 (both
  years) and P8 (2025)** — so their genetic test needs more mothers.
- All six model-pair populations (D1, D2, D3) and P6 and P42 of the
  genotyped sites are sampled well enough in both seasons; P41 (2025) and
  P40 (2026) at EO76 are not.
- At the three original genotyped sites (Part C, per EO) the mothers held
  would recover 96 % (EO67), 91 % (EO70) and 94 % (EO76) of the predicted
  alleles.
- Counting only seed cleaned so far, 20 populations of 2026 have no usable
  mother yet — the cleaning backlog, not a field gap.

<a id="fig-s4"></a>
![Figure S4](figures/step30_A_sampling_adequacy.png)

**Figure S4.** Is each population sampled well enough? One row per population, grouped by Bottleneck Lineage; open circle = 2025, filled = 2026; red = undersampled. **Panel A** — predicted share of the population's SRK alleles that its sampled mothers' seedlings would recover (bars = 95 % interval over simulated deme pools; dotted line = 90 %). **Panel B** — share of fertile plants in demes with no sampled mother (any value above 0 = undersampled). **Panel C** — mothers sampled (circle), mothers needed (|) and mothers with cleaned seed so far (×). Source: `step21_sampling_adequacy.py`. Data: [`step30_A_sampling_adequacy_2025.tsv`](tables/step30_A_sampling_adequacy_2025.tsv), [`_2026`](tables/step30_A_sampling_adequacy_2026.tsv); per deme: [`step30_A_sampling_adequacy_per_deme_2025.tsv`](tables/step30_A_sampling_adequacy_per_deme_2025.tsv).

---

## 4. Step 2 — selecting demes across the population (model populations)

Sampling every deme to its Step 1 target recovers 95 – 99.8 % of a
population's SRK alleles, but in a large population it is expensive (P34
needs 58 mothers in 2025 and 125 in 2026). Step 2 reaches the 90 % goal
with far fewer plants while spreading the sample across the population.

**Approach** (`step28_model_population_design.py`):

1. **Screen the demes.** A deme is unfit if it has a single fertile plant
   or if we hold fewer usable mothers than its Step 1 target (even after
   every sampled seed lot is cleaned).
2. **Select demes.** Each deme keeps its Step 1 target. Demes at least
   **100 m apart** — so selected demes are separate pollen pools and the
   sample covers the population's extent — are added one at a time, the
   one adding the most expected new SRK alleles per mother first, until
   the selection is expected to detect **90 % of the whole population's
   SRK alleles** (unfit demes included in the total).
3. **Guard against optimism.** Demes are chosen on one set of simulated
   allele pools; the stopping rule uses 20 independent sets, and the 5th
   percentile is reported. The selection is rerun on 20 fresh sets to
   measure how often each deme is chosen.
4. **Small populations.** When all fit demes together need at most 10
   mothers, all are taken (spacing is moot).
5. **Pick the mothers — at random, not by seed yield.** In each selected
   deme, a mother is **eligible** if her cleaned seed lot holds at least
   25 seeds, enough for the full 15 seedlings. One mother is drawn at
   random from the eligible mothers of every slick spot, and the
   remaining places are drawn at random from the deme's other eligible
   mothers. Lots with fewer than 25 seeds (germinated in full), then
   mothers whose seed is not cleaned yet, are used only when a slick spot
   or deme runs out of eligible mothers. The draw is reproducible (fixed
   random seed per population, year and deme).

**D1 plan (P36 small + P34 big, BL3).**

| Population | Year | Demes (unfit) | Selected demes | Mothers | Seeds → seedlings | SRK alleles detected (5th pct.) |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| P36 (small) | 2025 | 2 (0) | 2 | 4 | 100 → ~60 | 97 % (95 %) |
| | 2026 | 2 (0) | 2 | 4 | 100 → ~60 | 94 % (92 %) |
| P34 (big) | 2025 | 16 (4) | 7 | 40 | 1 000 → ~600 | 91 % (85 %) |
| | 2026 | 21 (2) | 7 | 45 | 1 110 → ~666 | 92 % (86 %) |
| **D1** | | | **18 deme-years** | **93** | **2 310 → ~1 386** | |

- **Unfit demes.** P34 2025: two single-plant demes, and demes 5 (158
  plants) and 12 (42) with too few mothers sampled. P34 2026: one
  single-plant deme, and deme 3 (50 plants).
- **Selected demes** (deme · mothers). P34 2025: 15 · 1, 8 · 3, 4 · 3,
  14 · 6, 6 · 8, 11 · 8, 2 · 11. P34 2026: 11 · 1, 6 · 3, 14 · 5,
  10 · 6, 18 · 9, 12 · 6, 2 · 15. P36: both demes in both years.
- **Stability.** 16 of the 18 deme-years are chosen in every rerun, the
  other two (P34 2025 deme 14, P34 2026 deme 12) in 55 – 60 % of reruns.
- **Why random rather than the best-seeded mothers.** Beyond 25 seeds,
  more seed does not give more seedlings, so seed count has no
  genotyping value — but ranking by it would pick the mothers that set
  the most seed, i.e. those that found compatible mates. That sample
  would overstate SRK diversity at the compatible end, understate mate
  limitation, and weaken the comparison of compatibility with seed set
  that D1 is built for (most of all in the small population). With the
  earlier best-seeded rule, 57 – 75 % of the selected mothers came from
  the top half of seed yield; with the random draw, 50 – 51 %, i.e. the
  selected mothers represent the full range of seed set among plants
  with enough seed to genotype.
- **Material.** 83 of the 93 mothers have cleaned seed ready to
  germinate; 10 (P34 2026, demes 10 and 12) need their seed cleaned
  first. One P34 2026 lot holds only 10 seeds, hence 2 310 seeds rather
  than 93 × 25.
- **Spacing.** A 200 m rule would cap P34 2025 at 82 % (only 6 demes that
  far apart); 100 m costs almost nothing compared with no rule (40 vs 43
  mothers in 2025, identical in 2026). P36's two demes are less than
  100 m apart, which is why the small-population rule keeps both.
- **The largest demes are not always chosen.** In P34 2026 the two
  largest demes (325 and 270 plants) are left out: smaller demes reach
  the 90 % goal more cheaply per mother. They can be forced in if the test
  needs deme-size contrasts within P34.

<a id="fig-s5"></a>
![Figure S5](figures/step30_D_model_design_curves.png)

**Figure S5.** D1 deme selection. One panel per population and season: the share of the population's SRK alleles expected to be detected (independent simulations) as demes are added in the order of Step 2, against the cumulative number of mothers; dotted line = 90 %; open square = every fit deme sampled. Source: `step28_model_population_design.py`. Maps of the selected demes (`figures/step30_D_model_design_demes_P34`, `_P36`) show slick spots at full resolution and are kept off the public repository; `./run_pipeline.sh design` writes them locally.

**Caveat.** The model treats demes as independent draws from the
species-wide SRK frequencies. Real demes have drifted further
(Project_outline § 5: 5 – 9 alleles observed where 11 – 32 were
predicted), so actual recovery will differ, and the spacing between
selected demes matters more in reality than in the model. D2 and D3 can be
planned the same way (`step28_model_population_design.py --pairs D2,D3`).

---

## 5. The locked recipe and the files to use

| Quantity | Value |
|---|---|
| Pollinator radius / drift unit | **50 m deme** (connected component of slick spots within 50 m) |
| Seeds to germinate per mother | **25** (= ceil(15 ÷ 0.60)) |
| Seedlings genotyped per mother | **15** |
| Step 1 — mothers per deme | Size-scaled: 90 % of the deme's SRK alleles detected, + ≥ 1 per slick spot |
| Unfit deme | Single plant, or fewer usable mothers held than its target |
| Step 2 — deme selection | Demes ≥ 100 m apart, best gain per mother, until 90 % of the population's SRK alleles are detected |
| Mothers within a selected deme | Random among cleaned lots with ≥ 25 seeds: one per slick spot, then the rest of the target; smaller lots, then uncleaned lots, only as fallback |
| Undersampled population | Left out for that season (< 90 % recovered, or a deme with no mother) |

**Files.**

- **Lab sheet for the model populations, keyed on occurrenceID** (germplasmIDs
  are still missing for the seed lots not cleaned yet):
  [`step30_D_model_design_lab_sheet.tsv`](tables/step30_D_model_design_lab_sheet.tsv)
  — one row per mother: occurrenceID, germplasmID, BL, population, EO,
  location, year, deme, slick spot, status ("germinate" / "clean seed
  first"), seeds held, seeds to germinate, expected seedlings and why the
  mother was selected.
- **Lab list for the model populations — the germplasmIDs to pull, deme by
  deme:**
  [`step30_D_model_design_germplasm_by_deme.csv`](tables/step30_D_model_design_germplasm_by_deme.csv)
  (one row per mother: germplasmID, occurrenceID, ready or "clean seed
  first", seeds held, seeds to germinate, expected seedlings, reason for
  selection). Deme screening and selection:
  [`step30_D_model_design_demes.tsv`](tables/step30_D_model_design_demes.tsv);
  totals: [`step30_D_model_design_summary.tsv`](tables/step30_D_model_design_summary.tsv).
- **All-population lab recipe (Step 1 targets, every deme):**
  [`step29c_partC_germplasmID_selection_2025.tsv`](tables/step29c_partC_germplasmID_selection_2025.tsv),
  [`_2026`](tables/step29c_partC_germplasmID_selection_2026.tsv).
- **Field-team recipe for a new season** (mothers per slick spot):
  [`step29c_sampling_frag_aware_per_event_2025.tsv`](tables/step29c_sampling_frag_aware_per_event_2025.tsv).
- **Sampling check per population:**
  [`step30_A_sampling_adequacy_2025.tsv`](tables/step30_A_sampling_adequacy_2025.tsv),
  [`_2026`](tables/step30_A_sampling_adequacy_2026.tsv).
- **Which slick spots share a pollen pool:**
  [`step29a_event_to_deme.tsv`](tables/step29a_event_to_deme.tsv).

**Sensitivity to the pollinator radius.** The 50 m deme is chosen on
pollinator biology. At 50 m the Step 1 design needs 695 mothers in 2025
(1 110 in 2026), 24 % (32 %) fewer than at 10 m; beyond 75 m it changes
little, and if the true radius were only 25 m the design would be 11 %
short (Project_outline Figure 1).
