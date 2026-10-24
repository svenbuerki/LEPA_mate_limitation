# Mate limitation versus SI escape in *Lepidium papilliferum* populations

**A compact overview for collaborators.** Framed on a single causal
chain (**fragmentation → genetic drift → mate limitation or SI
escape**); each step below is presented as Question → Approach → Result. This is
the project's main document. Earlier technical detail (formulas,
code, output filenames) is archived, no longer maintained, in
[`archive/Sampling_and_prediction.md`](archive/Sampling_and_prediction.md).

---

## Executive summary

The summary follows the order of the analyses: the breeding units
(§ 2), the genetic predictions (§ 3), the seed-yield test and how it
compares with the predictions (§ 4), the preliminary SRK genotypes
(§ 5), the sampling design that will test the predictions (§ 6), and
the proposed model populations (§ 7). § 8 discusses the results and
the next steps.

### 1. Overview

<span style="color:#777"><strong>Challenge.</strong></span> *Lepidium papilliferum*
(slickspot peppergrass, LEPA) is a federally listed tetraploid
Brassicaceae. This work covers the **Snake River Plain
populations**, in documented decline, as prioritised by current
stakeholder conservation needs; the Jarbidge Foothills populations
are out of scope and will be addressed separately. LEPA relies on
**sporophytic self-incompatibility (SI)**: a plant can only be
fertilised by a mate carrying different SRK alleles. Habitat
fragmentation shrinks the local breeding pool, and genetic drift
erodes the SRK alleles each plant needs to find a compatible mate.
Two reproductive failures can follow: **mate limitation** (too few
compatible mates, so seed set drops) or **SI escape** (the SI
barrier breaks down and plants self-fertilise, which keeps seed set
up but exposes the population to inbreeding).

<span style="color:#777"><strong>Goal.</strong></span> Establish a **population-level
restoration priority list**: which populations are mate-limited
and need intervention now, which show signs of SI escape, and
which are reproducing normally.

<span style="color:#777"><strong>Central hypothesis.</strong></span> **Fragmentation →
genetic drift → mate limitation**: drift in small, fragmented
demes erodes SRK alleles, fewer mates are compatible, and seed set
drops.

<span style="color:#777"><strong>Competing hypothesis — SI escape.</strong></span> A
population could instead escape mate limitation by losing SI and
self-fertilising. We test **how likely this competing hypothesis
is** with phenotypic data only — **seed yield and plant size**.

<span style="color:#777"><strong>How the two are told apart.</strong></span> **Seed
production reflects the reproductive status of a population.**
Plant size sets how many seeds a plant *could* make; whether it
makes them depends on how its pollination is working. Read against
a species-wide size → seed-yield expectation, **under-production
points to mate limitation**, production at expectation to intact
reproduction, and **over-production to SI escape**
(self-fertilisation adds seed beyond what outcrossing supplies).

<span style="color:#777"><strong>Approach.</strong></span> Two independent lines of
evidence on every population. **Genetic**: a drift model predicts
how many SRK alleles each deme holds and how often its plants are
compatible. **Behavioural**: observed seed yield is compared with
what plant size predicts; it needs no genotype, so it runs on all
populations now. SRK genotypes at six sites give a
first check of the genetic predictions, and seed–mother genotyping
in model populations will test them.

<span style="color:#777"><strong>Data.</strong></span> LEPA field database
(Genetic-Rescue-DB): 2025 + 2026 event and occurrence records, per-plant crown
and height, per-mother seed yield and germplasm records (**1 441
plants** with crown, height and seed yield). From the SRK
bioinformatics pipeline: the species-wide frequencies of **32
functional SRK alleles**, the empirical plant zygosity
distribution, and the genotypes of 246 plants (grown from wild seed or
sampled in the wild) at six sites (§ 5).

### 2. Breeding units — Bottleneck Lineages, populations and demes

<span style="color:#777"><strong>Methods.</strong></span>

- **Population** — connected component at 500 m on pooled 2025 +
  2026 event coordinates and historical location centroids. 500 m
  is ≈ 10 × the pollinator radius and well beyond LEPA's
  gravity-dominated seed dispersal, so populations are
  demographically independent. The population is stable across
  years, while the above-ground emergence of this annual plant
  flickers over a persistent seed bank (§ 8).
- **Bottleneck Lineage (BL)** — Ward's hierarchical clustering of
  population centroids into *k* = 5 lineages (fixed for continuity
  with the EO-level lineages; the silhouette optimum, *k* = 8, is
  only 0.03 higher). BLs are numbered by area, largest first
  (BL1 → BL5). Populations in a BL share a regional habitat-loss
  history, so BL controls shared drift history when comparing
  populations.
- **Deme** — 50 m connected component within one field season,
  nested in a population: the local pollen pool and the drift
  unit. The 50 m radius follows the pollinators' foraging range
  (halictid and small bees) and a Wright-neighbourhood argument
  (Wright 1943, 1946; Levin & Kerster 1974). A 10 – 200 m sweep
  shows how much the design depends on it ([Figure 1](#fig-1)).

<span style="color:#777"><strong>Results.</strong></span> The occurrences resolve into
**44 populations** ([Figure 2](#fig-2)) in **five BLs** (BL1 17
populations, BL2 9, BL3 10, BL4 5, BL5 3; [Figures 3](#fig-3) –
[5](#fig-5)). By 2025 → 2026 trajectory: **15 stable, 5 growth,
3 crash, 3 ambiguous, 7 above ground in 2025 only, 11 in 2026
only** ([Figures 6](#fig-6) and [7](#fig-7)). Demes are counted per
season, over the populations with plants above ground that year:

| Season | Populations above ground | Demes (mean per population) | Single-deme populations | Median deme size | Demes < 8 plants | Single-plant demes |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| 2025 | 33 | **99** (3.0) | 13 | 24 plants | 29 % | 6 % |
| 2026 | 37 | **189** (5.1) | 15 | 15 plants | 39 % | 11 % |

Below 8 plants a deme cannot carry one copy of every SRK allele in
the species pool; a single plant has no mate within 50 m
([Figure 8](#fig-8)). More than half of the populations breed as
several demes rather than one pool (2025: 20 of 33; 2026: 22 of 37).
The 2026-10-20 coordinate corrections removed four single-slick-spot
populations that existed only because of mistyped coordinates (48 → 44).

**Sampled vs recorded-only slick spots.** A few slick spots were
recorded but not sampled: the crew counted the fertile plants but
arrived after seed had dropped. These plants were part of that
season's mating pool, so they **count in the model** (census, demes,
predictions) but **supply no seed for experiments**:

| Population | Year | Fertile plants | In recorded-only slick spots | Demes with no seed |
|---|---|:-:|:-:|:-:|
| **P25** EO61 | 2026 | 523 | 176 (34 %) | 1 (20 plants) |
| **P4** EO30-3 | 2026 | 152 | 62 (41 %) | 1 (62 plants) |

### 3. Predicted SRK allele diversity and pollen compatibility

<span style="color:#777"><strong>Methods.</strong></span> For each deme, a
finite-population drift model draws the deme's SRK alleles (4 per
plant) from the species-wide frequencies. A population holds an
allele if any of its demes does ([Figure 9](#fig-9)). Pollen
compatibility is the share of random plant pairs within a deme that
can mate under **sporophytic Class I / Class II dominance with the
empirical zygosity**, averaged over demes weighted by size
([Figures 10](#fig-10) and [11](#fig-11)).

<span style="color:#777"><strong>Results.</strong></span>

- **SRK allele diversity** tracks deme size. The median population
  is predicted to hold 19 of the 32 species-wide alleles in 2025
  (21 in 2026), from 1 – 2 in a single-plant population to 29 – 32 in
  the largest. 13 populations in 2025 and 15 in 2026 are predicted to
  hold fewer than half (16). Each deme is simulated plant by plant
  with the observed zygosity — two thirds of LEPA plants carry a
  single SRK allele in all four copies — which lowers these numbers
  against an earlier version that drew each plant's four copies
  independently (median 26 – 27).
- **Pollen compatibility barely varies**: median 0.78 in both years
  (species mean 0.78). Only populations of 1 – 3 plants fall into the
  struggling band (0.26 – 0.52) or below: P31 (EO24-1), P32 (EO24) and
  P33 (EO24-2, one plant, 0.09) in 2025 and P33 in 2026, with P30
  (EO24-7) and P38 (EO52) at its edge (0.52). Three mechanisms buffer compatibility
  against allele loss — Class I dominance (26 of 32 alleles), the
  species' current class imbalance, and tetraploid homozygosity — so
  the model predicts a shortfall only in the smallest demes.

### 4. Seed yield by plant size — and how it compares with the predictions

<span style="color:#777"><strong>Methods.</strong></span>

1. **Species-wide curve.** Cross-validation picks crown size as the
   best predictor of seed yield (log scale; R² = 0.38). The curve is
   fitted on all plants from both years pooled
   ([Figure 12a](#fig-12a)).
2. **Per population and year.** Each plant's residual (observed −
   expected) is tested against 0 (one-sample t-test, corrected for
   multiple testing with Benjamini–Hochberg). The fold of expectation
   is observed / expected seed yield ([Figure 12b](#fig-12b)).
3. **Year check.** Within populations, 2026 plants made 15 % fewer
   seeds for their size than 2025 plants (fold 0.85, p = 0.005). The
   test is therefore run twice: with the pooled curve, and with this
   year effect as an offset on 2026. Only statuses that agree are
   ranked; the others are flagged **year-sensitive** (they depend on
   whether the 2026 drop is a year effect or a real decline).
4. **Population status.** The two years combine into one status:
   *mate-limited* (significant under-production in ≥ 1 year, never
   over), *switching* (under in one year, over in the other),
   *candidate mate-limited* (fold < 0.8, not significant),
   *SI-escape candidate* (significant over-production, never under),
   *reproduction intact*, or *insufficient data* (no year with ≥ 3
   plants).
5. **Comparison with the predictions.** Rank correlation (Spearman)
   between each population-year's fold of expectation and its
   predicted pollen compatibility and SRK diversity
   ([Figure 12c](#fig-12c)).

<span style="color:#777"><strong>Results — the restoration priority list</strong></span>
([`step30i_conservation_priority.tsv`](tables/step30i_conservation_priority.tsv)):

| Rank | Status | n | Populations |
|:-:|---|:-:|---|
| 1 | **Mate-limited** | **6** | P1 (EO112, new to our sampling in 2026), P4 (EO30-3), P8, P10 (EO27 cluster) (BL1) · P32 (EO24, BL3) · P42 (EO70, BL5) |
| 2 | **Switching** (over in 2025 → under in 2026) | **1** | P19 (EO8, BL2) |
| 3 | Candidate mate-limited | 3 | P3 (EO30-4), P6 (EO67) (BL1) · P20 (EO26-2, BL2) |
| 4 | **SI-escape candidate** | **7** | P2 (BL1) · P18, P23, P24, P25 (BL2) · P37, P40 (BL4) |
| 5 | Reproduction intact | 8 | — |
| — | **Year-sensitive** | **9** | mate-limited with the pooled curve: P16 (EO27-1 + EO27RT; switching with the offset), P36 (EO25-B; candidate) · SI-escape with the pooled curve, intact with the offset: P26, P29, P34, P35, P43 · intact → candidate: P13, P27 |
| 6 | Insufficient data | 3 | P7, P31, P33 |

The most severe shortfall is P10 EO27-3 (fold 0.44 in 2025, 0.02 in
2026); the strongest over-production is in BL2 (P18 2.26, P23 2.19)
and at P37 (2.47 in 2026) and P40 (2.22).

<span style="color:#777"><strong>Results — seed yield vs the predictions.</strong></span>
**Seed yield does not track predicted pollen compatibility** (ρ = 0.17,
p = 0.27, 45 population-years; 2026 alone ρ = 0.45, p = 0.08; 2025
ρ = 0.09, n.s.). Predicted SRK diversity and census do not correlate
either (ρ = 0.17). An earlier version of the drift model, which drew
each plant's four allele copies independently, gave ρ = 0.38
(p = 0.01); simulating plants with the observed zygosity removed that
signal. **The mate-limited populations sit at the species mean**
(0.72 – 0.79; P8, P10, P42 at 0.77 – 0.79): the drift model, which
draws each deme once from the species-wide frequencies, does not
predict where mate limitation occurs. Seed yield is the
direct read-out of reproductive status, and the populations where
the two disagree are the first targets for genotyping (§ 5, § 8).

### 5. Preliminary results — the genotyped populations

SRK genotypes exist for plants grown from wild seed (and some wild
leaf tissue) at six EOs. Each plant is placed in a population through
its seed lot (wild mother → her slick spot), its sample location or
its tissue position. That gives **six sites, 246 plants** (§ C.0.a,
[Figure 13](#fig-13)): three EO sites — **EO67 = P6, EO70 = P42,
EO76 = P40 + P41** (compared as one site) — and three populations of
the split EOs with ≥ 10 plants — **P34 (EO18), P13 and P16 (EO27)**.
EO25 cannot be split between its two populations and stays out.

- **SRK allele diversity is far below the prediction at five of the
  six sites**: 5 – 9 alleles observed against 22 – 29 predicted; only
  the small P6 matches (7 observed, 6.5 predicted). A radius sweep at
  EO70 (P42) and EO76 (P40 + P41) ([Figure 15](#fig-15)) shows the gap stays flat down to
  10 m, so it is not a deme that is too wide: **local demes have
  drifted further than one draw from the species-wide frequencies
  assumes**.
- **Pollen compatibility follows the prediction at P6** (0.70
  observed, 0.62 predicted, within its wide interval 0.39 – 0.74)
  **and EO76** (0.72 / 0.78), **falls below it at P42** (0.53 vs
  0.78), **P34** (0.67 vs 0.79) **and P13** (0.44 vs 0.79), and is
  **above it at P16** (0.84 vs 0.79; 0.77 with P16's own, unusually
  multi-allele zygosity). The decomposition ([Figure 14](#fig-14))
  finds three different causes: one dominant allele within a class at
  P42 (FG024), the loss of almost all Class I alleles at P13 (4 % of
  its SRK copies), and genotype zygosity at P16. **Which alleles are
  lost and how they combine in plants, not the raw allele count,
  drives compatibility.** P13 and P16 rest on 12 and 15 plants (7
  mothers each) and are indicative.
- **Agreement with seed yield**: P42 is mate-limited in both years
  (robust to the year offset) — a compatibility gap and a seed
  shortfall in the same population. P40 over-produced in 2025 with
  only 9 of 32 alleles before it dropped to 23 plants in 2026,
  consistent with selfing in a pool short of compatible mates.
  Genotyped seedlings are counted as independent plants; siblings of
  one mother are not yet corrected for.

### 6. Sampling design

Testing the predictions needs seed–mother genotyping (each seedling
carries 2 paternal SRK alleles from its mother's pollen pool). The
sampling design and the check that our material can test the
predictions are in a dedicated document, [`Sampling_design.md`](Sampling_design.md). In short, a
**two-step plan recovers 90 % of each population's predicted SRK
diversity**: (1) within each deme, enough mothers (25 seeds → 15
seedlings each) to detect 90 % of the deme's SRK alleles, sized to the
deme, with at least one mother per slick spot; (2) across the
population, demes at least 100 m apart are added until 90 % of the
population's SRK alleles are recovered with the fewest mothers. With
the material we hold, 27 of 33 populations (2025) and 32 of 37 (2026)
are sampled well enough.

### 7. Proposed model populations

We propose a **within-BL small-vs-big** pair: a small and a big
population in the same Bottleneck Lineage, both present in both
years with at least 20 fertile plants summed across years. Holding
BL constant means a difference between the two is a small-vs-big
signal, not shared drift history. BL1, BL2 and BL3 each hold a
candidate pair; BL4 and BL5 do not. All six populations are sampled
well enough in both seasons (§ 6).

| Option | BL | Small — fertile plants (2025 + 2026), mothers sampled (2025 + 2026) | Big — same | Size ratio | Mothers needed, every deme (2025 + 2026) | Comment |
|:---|:---:|:---|:---|:---:|:---:|:---|
| **D1 ★** | BL3 | **P36 EO25-B** · 26 plants · **7 + 8** mothers | **P34 EO18-7 + EO18-8** · 1 673 plants · **115 + 239** mothers | 64× | **4 + 3 / 59 + 98** | Largest mother pool; both sides sampled well enough in both years. |
| D2 | BL2 | P21 EO26-3 + EO26 · 191 plants · 18 + 52 mothers | P19 EO8 · 876 plants · 37 + 178 mothers | 5× | 9 + 28 / 22 + 90 | Neither side is stable (P21 ambiguous, P19 growth); smallest size contrast. |
| D3 | BL1 | P6 EO67 · 22 plants · 4 + 6 mothers | P16 EO27-1 + EO27RT · 1 991 plants · 55 + 251 mothers | 90× | 2 + 2 / 53 + 136 | Biggest size contrast and both sides stable, but the small side has only 10 mothers. |

**D1 is recommended** (★). It has a sharp 64× size contrast, both
sides in BL3, and both sides sampled well enough in both years
([`Sampling_design.md`](Sampling_design.md) § 3). With the pooled curve
**P36 is mate-limited** (fold 0.52 in 2026) and **P34 over-produced
in 2025** (1.28; 0.91 in 2026, not significant), the contrast D1 is
built on. Both calls are **year-sensitive**: with the 2026 year offset
P36 becomes a candidate and P34 intact, so the pair's seed-yield
contrast is not yet firm — seed–mother genotyping will test it
directly. **Caveat:** P34's census is not constant (523 → 1 150
plants, 16 → 21 demes); no large population has both a constant
census and seed yield at expectation in both years. D2's big side, P19, is the one
switching population and the natural follow-up.

**What D1 is for.** D1's first job is **SRK diversity**: measuring
how many alleles a small and a big population of the same lineage
actually hold, deme by deme, to calibrate the drift model (§ 8). With
91 mothers it is too small to test mate limitation mother by mother on
its own; that test needs mothers from many populations
([`Sampling_design.md`](Sampling_design.md) § 1).

**D1 sampling plan.** Selecting spatially separated demes (Step 2 of
the sampling design) recovers about 90 % of each D1 population's SRK
alleles with **91 mothers over the two years** (2 260 seeds to
germinate, ≈ 1 360 seedlings), all with cleaned seed. In 2025 P34's
largest deme holds fewer mothers than its target and is used with the 9
it holds, the one adjustment to the plan; the selected demes and the
germplasmIDs to pull are in [`Sampling_design.md`](Sampling_design.md) § 4.

Full `populationID → locationCode` map in
[step30g_populations_classified.tsv](tables/step30g_populations_classified.tsv).

### 8. Discussion and next steps

**Mate limitation is concentrated, and it is not where the drift
model expects it.** The EO27 cluster (P8, P10; P16 year-sensitive),
P42 EO70, P32 EO24 and two small late-season BL1 populations (P4, P1)
under-produce seed, yet most are large and predicted at the
species-mean compatibility; across all populations seed yield does
not follow predicted compatibility (ρ = 0.17, n.s.). The
genotyped sites point to the reason: local demes have lost far more
SRK alleles than one draw from the species-wide frequencies
predicts, and at P42, P34 and P13 that loss has skewed the allele
frequencies or the class balance enough to cut compatibility. Drift has acted longer, or on smaller
effective sizes, than our model assumes. Not every shortfall need be
genetic, though — pollinators, resources, inbreeding depression and
flowering time can also depress seed set. Seed–mother genotyping
separates the two: low observed compatibility points to drift,
normal compatibility to another cause.

**SI escape is a lineage-scale signal.** Five BL2 populations
over-produce for their size in 2025 (P18, P23, P24, P25 robustly, P26
year-sensitive) plus the switching P19. BL2 is the mid-area, high-connectivity lineage, so this is not
a small-population artefact. If confirmed, SI escape is a
conservation concern in its own right: selfing erodes the SI system
and exposes populations to inbreeding depression.

**The seed bank buffers drift and makes populations swing.** LEPA
seeds persist in the soil for up to 12 years: each year about 5.8 %
of a seed cohort germinates and 2.6 % dies (Meyer, Quinney & Weaver
2006, *Am. J. Bot.* 93: 891 – 902). The plants above ground each year
are therefore drawn from about 12 seed cohorts, and the effective
generation time is about 6.5 years. Two consequences follow:

- **Across years the bank slows allele loss** — a population keeps
  alleles that are absent above ground in any one year (the storage
  effect). The low diversity of the genotyped plants must then come
  from past drift in the local bank, not from one bad year.
- **Within a year the mating pool is a small, changing draw** from
  the bank, so census, compatibility and seed set can swing between
  years: P19 (fold 1.74 → 0.60) and P34 (1.28 → 0.91). Because only
  ≈ 6 % of the 2025 seed germinates in 2026, selfed 2025 seed cannot
  explain most of a 2026 drop; turnover of the emerging plants, or
  weather, are the likelier causes.

**Limitations.** Seed counts flag candidates; confirming the
mechanism (outcrossing vs selfing) needs seed–mother genotyping,
since selfed seedlings carry only maternal SRK alleles. 2026 seed set
is 15 % below 2025 within populations; with two years we cannot tell a
year effect from a real decline, so nine statuses are year-sensitive.
2026 calls should be re-run once the germplasm clean-up is complete. The 50 m
deme is a geographic upper bound on realised gene flow. The
seed-bank rates come from one study's sites and are averages; they
vary with weather.

<span style="color:#777"><strong>Next steps.</strong></span>

1. Re-run the seed-yield test and the sampling check once the 2026
   germplasm records are complete, and add 2027 to separate the year
   effect from population decline.
2. Add the seed bank to the drift model: a 12-cohort simulation per
   deme that predicts allele loss with and without the bank, the
   expected similarity of a population's SRK frequencies between
   years, and whether seed-bank turnover alone can produce swings
   like P19's.
3. **Calibrate the drift model.** Real populations hold far fewer SRK
   alleles than draws from the species-wide frequencies predict (§ 5),
   and predicted compatibility does not track seed yield (§ 4). Fit how
   far demes have drifted — an effective size well below the census,
   or allele frequencies per Bottleneck Lineage instead of species-wide
   — to the genotyped sites, then re-predict every population.
4. Genotype seedlings and their mothers in the model populations
   (§ 7): their SRK diversity, measured in a small and a big
   population of the same lineage, is the data the calibration needs;
   they also measure compatibility and the selfing rate directly, and
   2025 vs 2026 lots test the seed-bank predictions.
5. Collect more mothers in the undersampled populations, starting
   with the mate-limited P10 and P8.
6. Extend seed–mother genotyping to the SI-escape candidates,
   starting with BL2. Testing mate limitation mother by mother needs
   many mothers: D1's 91 detect a strong effect (half the compatible
   mates → 50 % less seed) only 38 % of the time, while the ~1 060
   mothers held across every population sampled well enough detect a
   30 % effect 79 % of the time ([`Sampling_design.md`](Sampling_design.md) § 1).

---

## The central hypothesis — a causal chain

*Lepidium papilliferum* (slickspot peppergrass, LEPA) is a tetraploid
Brassicaceae with **sporophytic self-incompatibility (SI)** — a plant
rejects pollen carrying any of the SRK alleles the pollen parent
expresses on its own stigma. In small, patchy populations we
hypothesise a single causal chain that drives reproductive failure:

> **Fragmentation → genetic drift → mate limitation.**
>
> Habitat fragmentation shrinks the effective deme inside each
> population (fewer plants within ~50 m pollinator flight of each
> other). Small effective demes intensify genetic drift on
> SRK, which erodes local SRK diversity and skews local Fg
> composition. The eroded and skewed local pool reduces the fraction
> of pollen a mother is compatible with — **mate limitation** — which
> reduces per-mother seed set.

The pipeline evaluates the chain in order — **fragmentation first**,
then its drift consequences, then its mate-limitation consequences —
because fragmentation is the physical driver upstream of everything
else. It is where population size enters the model (via the
per-deme `component_N_fertile`); a census of 1 000 plants split
across 20 disconnected slick spots behaves like 20 small drift-prone
demes of ~50 plants each, not one pool of 1 000.

| Link | What we measure | Steps that produce it |
|---|---|---|
| **1. Fragmentation** | 50 m deme structure per population (count + size per deme) | Step 29b, Step 29c, Step 29d (§ A.5, § B.2 of long doc) |
| **2. Genetic drift on SRK** | Predicted local allele pool size + frequency composition **per deme**, aggregated to the population | Step 30 Phase A (§ A.6) |
| **3. Mate limitation** | Predicted per-population random-mating pollen compatibility `P_compat` under sporophytic Class I / II SI | Step 30 Phase A (§ A.7, § A.8); Step 30c empirical validation on adult SRK genotypes (§ C.0.a) |
| **4a. Behavioural seed-set check** (**no SRK genotype needed — runs on all 44 populations**) | Observed seed yield per plant vs the species-wide size → yield allometry; fold < 1 at FDR < 0.05 = under-producing (mate limitation); fold > 1 at FDR < 0.05 = over-producing (SI escape); fold ≈ 1 = reproduction intact | Step 30i (§ 30.3) |
| **4b. Mate-limitation regression** (full test — needs seed genotypes) | Observed per-mother seed set regressed on **observed** `P_compat` from seed-father genotypes, decomposed into drift (β₁) and fragmentation (β₂) channels | Step 30 Phase B (§ C.1 of long doc) |

**SI escape.** Where SI has broken down, plants self-fertilise and
set more seed than their size predicts. The seed-yield test (link
4a) flags these populations as SI-escape candidates; seed–mother
genotyping confirms selfing. Individual self-compatible (SC / pSI)
plants are called separately in the Canu-amplicon SRK genotyping
(Phase 4 Step 22b) and enter the framework as prior information
via the empirical zygosity distribution.

---

## Scope and data

- **Species.** Tetraploid (2n = 4x). Every plant carries 4 SRK
  allele copies. Sporophytic Class I / Class II SI: Class I strictly
  dominant within a plant; between-class crosses always compatible;
  within-class crosses compatible only when expressed sets are disjoint.
- **Zygosity.** Canu-amplicon Step 23 shows ~66 % of adult plants
  carry a single distinct functional SRK identity (homozygote-like),
  32 % carry two, ~2 % carry three. This empirical distribution is
  the mother- and father-drawing prior in every P_compat calculation.
- **Data source.** `LEPA_SQL.db` — wild in-situ occurrences,
  filtered to records with coordinates inside LEPA's Idaho range.
  The 2025 single-year scope used by the Phase 4 adult SRK
  validation (§ C.0 and the figures below): **237
  events in 33 populations (99 demes), 765 mothers with seed
  records.** The 2025 + 2026 pooled scope used by the
  Population framework (§ Populations) and the across-year dynamic
  (§ Across-year population dynamic): **44 populations, 559
  slick spots after 10 m cross-year matching, 288 per-year demes.**
- **Species-wide SRK prior (P1).** The **32 functional SRK allele
  groups (Fgs)** identified across LEPA, together with their
  species-wide empirical frequencies — a 32-slot probability
  distribution built from the Canu-amplicon L1 carrier inventory
  (49 sequence alleles collapsed into 32 Fgs across 263 individuals).
  "Draw an allele from P1" means pick one of the 32 identities with
  probability equal to its share of the inventory. Dominant Fg =
  FG001 (**41 %** of P1 mass). **Class I = 26 Fgs (~35 % of P1
  mass), Class II = 6 Fgs (~65 %)** — this inverts the "common =
  Class I" Brassica frequency pattern, which assumes neutrality. In
  LEPA genetic drift elevates FG001-006 regardless of class. See
  § A.7.2 of the long doc for the drift note; a phylogenetic
  reassignment is tracked as future work.

---

## Key concepts and terminology

Every per-population quantity in this doc is derived from a nested
spatial hierarchy. The top three levels are **framework constructs**
built on top of the LEPA DB; the two DB-level entries below them are
**quoted verbatim from the LEPA DB `Terms` table** (the canonical
glossary that ships with `LEPA_SQL.db`). Reading from the biggest
unit down to the individual plant:

- **Species** — *Lepidium papilliferum* on the Snake River Plain,
  summarised at the SRK level by the 32-Fg species-wide prior P1.
- **Bottleneck Lineage** (BL) — a Ward's hierarchical cluster of
  populations sharing an evolutionary context. 5 BLs in LEPA
  (*k* = 5 fixed, on the silhouette plateau, § Bottleneck Lineages); ordered BL1 … BL5 by
  convex-hull area DESC → within-BL connectivity DESC. A
  derived concept, not a DB term.
  - **Population** (`populationID`) — a 500 m-separated connected
    component on the pooled 2025 + 2026 event coordinates +
    historical locationID centroids. 44 populations in the current
    data. 500 m is ≈ 10 × the pollinator radius and well beyond
    LEPA's seed dispersal, so between-population gene flow is ≈ 0
    on an ecological timescale. The population is the **stable
    spatial patch** that persists even in years when above-ground
    emergence is zero at a given slick spot. A derived concept,
    not a DB term; carries a comma-separated legacy label
    `{locationCode}_{locationID}` for human readability (e.g.
    `EO27-1_11, EO27RT_12` for the merged population P16).
    - **Deme** (`component_id_50m`) — a group of events whose
      plants sit within 50 m of each other **within one year**,
      nested inside a population. **Plants in the same deme share
      a pollen pool; plants in different demes — even in the same
      population — do not.** This is the within-year drift unit on
      which every per-population prediction is built. A
      derived concept, not a DB term.
      - **Event** (`eventID` / `occurrenceID`) — DB `Events` table:
        "**an 'Event' refers to an occupied slick spot within a
        Location**" (Darwin Core `dwc:eventID`). Each event has its
        own census of fertile plants and its own coordinates. EventIDs
        are NOT stable across years — field crews assign fresh
        barcodes each season, so the framework matches slick spots across
        years with a 10 m haversine buffer before the population
        graph is built.
        - **Mother plant** (`germplasmID`) — an individual plant
          already collected and stored in the LEPA DB, sitting at
          one specific event.

**Legacy DB label — Location** (`locationID` / `locationCode`) —
DB `Locations` table: `locationID` = "Location Unique Barcode #"
(Darwin Core `dwc:locationID`); `locationCode` = the EO code
("Report the unique EO # where the sampling is conducted
(e.g., EO38)"). A single `locationID` can appear inside one
population (1:1 mapping — the common case) or get merged with
sibling locationIDs when ≤ 500 m apart (e.g. P16 = EO27-1_11 +
EO27RT_12). The refinement that split EOs ≥ 500 m apart
into dash-suffixed `locationCode`s (EO24, EO24-1, EO24-2, EO24-7;
EO27, EO27-1, EO27-3, EO27RT) is preserved in the crosswalk
(`step29a_population_crosswalk.tsv`, local only: it holds event coordinates)
so every legacy reference stays traceable through the new
population frame.

Everything else in the pipeline is a count or a derived number
sitting on top of this hierarchy:

| Concept (full English name) | Code identifier | Rooted at | Definition |
|---|---|---|---|
| Fertile plant census | `total_n_fertile` | population | All fertile plants in a population, summed across every event in a given year. The biological potential, with no within-year spatial filtering. |
| **Deme size** | `component_N_fertile` | deme | The fertile plants that share a single 50 m pollen pool in one year — **this is the drift unit for both diversity and pollen-compatibility prediction**. Each 50 m connected component within a population × year = one deme. A population can hold one deme (fully connected in a given year) or several (fragmented); see § Delineating the deme. |
| Fragmentation-aware mother target | `M_frag` | event, derived from demes | For each event, the number of mothers to sample so that each 50 m deme's seedlings detect 90 % of its SRK alleles on average (species-wide frequencies, drift), with a ≥ 1-per-event maternal-genotype floor. Sums across events to the population-level `M_frag_aware`. |
| Tetraploid per-mother seed cap | 15 seeds/mother | mother plant | Each seed contributes 2 paternal allele draws from the local pollen pool. 15 seeds/mother is the per-mother floor that gives a 90 % chance of seeing every allele in her deme's pollen pool. |
| Species prior | `P1` | species-wide | The **32 Fgs identified across LEPA plus their empirical species-wide frequencies** — a 32-slot probability vector that sums to 1, built from the Canu-amplicon L1 carrier inventory. Common Fgs (e.g. FG001 at 41 %) have a large slot; rare ones have a small slot. Every per-population prediction draws alleles from P1, so small populations lose the rare Fgs to drift by chance. |

**Why demes matter in one sentence.** Every per-population
prediction in this doc — SRK diversity, pollen compatibility,
mother allocation — is built **deme-by-deme**, because a 50 m
connected component is what trades pollen within one year. Each
deme gets its own `component_N_fertile`; population-level numbers
are the **set union** across demes for diversity (Fgs are a set)
and the **size-weighted mean** across demes for pollen
compatibility (a continuous rate). A population with 500 fertile
plants spread across 20 isolated slick spots behaves like 20 small
drift-prone demes, not one pool of 500.

### How the five quantities chain together

1. **Population → events.** Raw census per year `total_n_fertile`
   = Σ event `n_fertile` across the population's events in that
   year (step28_events_spatial_neighborhood.tsv).
2. **Events → 50 m demes.** Within each population × year, connect
   events whose fertile plants sit within ≤ 50 m; the connected
   components are the demes
   (step29a_event_to_deme.tsv).
3. **Demes → effective deme size.**
   `component_N_fertile_c` = Σ_{events ∈ c} `n_fertile_e` — the
   plants that actually share one 50 m pollen pool and the drift
   unit for both predictions.
4. **Effective deme → SRK diversity.** For each deme,
   draw `4 × component_N_fertile_c` alleles from P1 → the deme's
   present Fgs. Population pool size = |union of per-deme Fg sets|.
   How well our sampling recovers this pool is covered in
   [`Sampling_design.md`](Sampling_design.md).
5. **Effective deme → pollen compatibility.** For each deme,
   sample mothers from its local frequencies under the empirical
   LEPA zygosity and compute sporophytic Class I / II P_compat per
   mother. Population P_compat = size-weighted mean of per-deme
   P_compat.

The whole causal chain **fragmentation → drift → mate limitation**
enters at step 3 and comes out at steps 4 and 5: fragmented
populations have small demes → small drift-prone pools → fewer
Fgs and lower pollen compatibility than a one-pool population of
the same raw census would predict.

---

## The question this doc builds toward

**How many seeds per mother should the lab genotype for Part C
testing?** That single number is the practical output of the whole
framework, and it depends on three things — none of which can be
skipped:

- **Fragmentation** — how many adults actually share a pollen
  environment inside each population.
- **SRK diversity** — how many Fg alleles the population holds
  locally under fragmentation-driven drift.
- **Pollen compatibility (P_compat)** — how much of the local pool
  a mother is compatible with under sporophytic Class I / II SI.

The next sections build those three quantities up in order:
Populations (§ Populations) define the stable spatial frame,
Bottleneck Lineages (§ Bottleneck Lineages) stratify them by
shared drift history, the across-year dynamic (§ Across-year
population dynamic) classifies each population's trajectory, and
only then do we delineate demes inside each population and run the
Phase A predictions. The sampling design — mother count per
population, seed count per mother — is derived at the end, once
the causal chain is on the table.

---

## Populations — the stable spatial frame (Step 29a)

**Why a population level above the deme.** *L. papilliferum* is an
**annual** whose above-ground presence at a slick spot can switch
on and off between years depending on seed-bank germination and
growing-season conditions. A single-year "location" footprint is
therefore unstable — a slick spot can show plants in 2025 and not
in 2026, or vice versa, while the **underlying seed-bank patch
persists**. The framework therefore separates the stable **population**
(the patch, demographically independent spatial unit) from the
within-year **deme** (the drift unit inside it), following
classical metapopulation theory (Levins 1969; Hanski 1998;
Freckleton & Watkinson 2002).

| Level | Operational definition | Biological meaning |
|---|---|---|
| **Population** | 500 m-separated connected components, built from pooled **2025 + 2026 event coordinates + historical locationID centroids** | Demographically independent spatial patch; between-population pollen (≈ 10 × the 50 m flight radius) and seeds (>> LEPA's gravity-dominated dispersal) do not flow on ecological timescale |
| **Deme** | 50 m connected component **within one year's event set**, nested inside a population | Operational deme / Wright-neighbourhood scale — the within-year pollen pool and drift unit (see below) |
| **Event** | one occupied slick spot in one year | Darwin Core `dwc:eventID` — the raw observation record |

**Why 500 m.** (i) **≈ 10 × the 50 m pollinator radius** — far
outside the ecologically-realised pollen neighbourhood. (ii) **>>
LEPA seed dispersal** — Brassicaceae in low-stature arid
vegetation disperse by gravity + short-distance wind on the order
of metres. (iii) **Consistent with the existing within-EO 500 m
rule** already used to split DB locationCodes 
(5 EOs → 16 locationCodes under the current rule) —
promoting it to the formal population threshold removes a loose
end rather than adding a new parameter.

**Operational algorithm.** Field crews assign fresh event IDs each
season, so first we match slick spots across years with a **10 m
coordinate buffer** (absorbs typical GPS drift); then pool
slick spot centroids + historical locationID centroids and build a
500 m haversine graph whose connected components = populations.
Each DB locationID, locationCode, slick spot, and
(year, eventID) is preserved in a crosswalk TSV so legacy
references stay traceable. Within each population × year, rebuild
the 50 m deme partition — the within-year drift unit used by every
downstream prediction.

**Downstream.** The old "location" becomes a legacy label.
SRK diversity and pollen compatibility become **year-resolved per
population**: within one population, 2025 and 2026 above-ground
samples are two independent draws from the same seed bank, so
systematic drift between them directly estimates per-population
drift beyond the species-wide prior — a stronger test of H2
(§ C.0.c) than the current one-year sweep, built on this project's
own data. The step also outputs a **seed-cleaning priority
ranking** that flags one large and one small candidate population
(≥ 2 re-visited slick spots in both years, consistent above-ground
trend) as the preliminary across-year analysis targets — the
2026 seed-cleaning queue is sized around that priority list.

### First run (2025 + 2026 pooled)

Site visitation effort was equal across years; the extra 2026
occurrence records (mother plants collected for seed banking) are
not used at this stage, so raw `N_fertile` sums are the appropriate
across-year axis.

- **44 populations** across the Snake River Plain, built from
  559 slick spots (from 723 raw events after 10 m cross-year
  matching) + 52 historical locationID centroids.
- **26 populations present in both 2025 and 2026** / 7 only in
  2025 / 11 only in 2026. 288 per-year demes (99 in 2025, 189 in
  2026).
- **55 of 559 slick spots (9.8 %) were re-visited in both years** —
  the raw signal of the annual's year-to-year above-ground
  flickering.

**Legacy locationCode mapping.** 28 locationCodes map 1:1 to a
population; 9 span several populations; 4 populations aggregate ≥ 2
legacy locationCodes (P16 merges EO27-1 + EO27RT, the clearest case).

<a id="fig-2"></a>
![Figure 2](figures/step29a_populations_overview.png)

**Figure 2.** The 44 populations across the Snake River Plain — 2025 + 2026 above-ground occupancy. One dot per population centroid, size ∝ total `n_fertile` across both years, coloured by occupancy pattern: purple = present in both years (26 of 44 populations), blue = 2025 only (7), orange = 2026 only (15). 44 populations total built from 562 slick spots (10 m cross-year matching) plus historical locationID centroids; 500 m connectivity graph on the pooled point set. 9.3 % of slick spots were re-visited at the 10 m scale across both years — the raw signal of the annual's year-to-year above-ground flickering. Source: `step08_population_maps.py`.

The two downstream sections pick up from here: § **Bottleneck
Lineages** below groups the 44 populations into the five BLs
that structure every per-population figure, and § **Across-year
population dynamic** classifies each population's 2025 → 2026
trajectory and names the first candidate test pairs.

### populationID ↔ locationCode(s) lookup

The inline table below is the authoritative populationID ↔
locationCode(s) + BL + occupancy + trend-class map for the 48
populations. Every per-population figure that uses
short `P{N}` row labels resolves here; the same information is
also in the TSV
[`step30g_populations_classified.tsv`](tables/step30g_populations_classified.tsv)
(which carries additional columns — N_fertile per year,
predicted diversity, predicted pollen compatibility, cross-year
flags). Rows ordered by BL (BL1 → BL5) then populationID.

| populationID | BL | locationCode(s) + locationID(s) | occupancy | trend class |
|:-:|:-:|---|:-:|:-:|
| P1 | BL1 | EO112_53 | 2026 only | 2026_only |
| P2 | BL1 | EO30-1_13 | both years | stable |
| P3 | BL1 | EO30-4_52 | 2026 only | 2026_only |
| P4 | BL1 | EO30-3_51 | 2026 only | 2026_only |
| P5 | BL1 | EO30-2_42 | 2026 only | 2026_only |
| P6 | BL1 | EO67_39 | both years | **stable (SMALL)** |
| P7 | BL1 | EO72-2_40 | 2025 only | 2025_only |
| P8 | BL1 | EO27_9 | both years | stable |
| P9 | BL1 | EO27_9 | 2026 only | 2026_only |
| P10 | BL1 | EO27-3_10 | both years | crash |
| P11 | BL1 | EO27_44 | 2026 only | 2026_only |
| P12 | BL1 | EO27-5_43 | 2026 only | 2026_only |
| P13 | BL1 | EO27-1_37 | both years | stable |
| P14 | BL1 | EO27-1_45 | 2026 only | 2026_only |
| P15 | BL1 | EO27-1_46 | 2026 only | 2026_only |
| P16 | BL1 | EO27-1_11, EO27RT_12 | both years | **stable (BIG)** |
| P17 | BL1 | EO27_50 | 2026 only | 2026_only |
| P18 | BL2 | EO8_29, EO8_47 | both years | growth |
| P19 | BL2 | EO8_27, EO8_28 | both years | **growth (BIG)** |
| P20 | BL2 | EO26-2_35 | 2025 only | 2025_only |
| P21 | BL2 | EO26-3_32, EO26-3_33, EO26-3_34, EO26_48 | both years | **ambiguous (SMALL)** |
| P22 | BL2 | EO26-4_36, EO26_49 | both years | stable |
| P23 | BL2 | EO26-3_31, EO26-3_32 | both years | stable |
| P24 | BL2 | EO26-1_30 | 2025 only | 2025_only |
| P25 | BL2 | EO61_38 | both years | stable |
| P26 | BL2 | EO29_8 | both years | stable |
| P27 | BL3 | EO48_7 | 2025 only | 2025_only |
| P28 | BL3 | EO32_6 | both years | growth |
| P29 | BL3 | EO32_6 | both years | stable |
| P30 | BL3 | EO24-7_25 | 2025 only | 2025_only |
| P31 | BL3 | EO24-1_22 | 2025 only | 2025_only |
| P32 | BL3 | EO24_24 | both years | growth |
| P33 | BL3 | EO24-2_23 | both years | ambiguous |
| P34 | BL3 | EO18-7_15, EO18-7_16, EO18-7_17, EO18-7_19, EO18-8_18 | both years | **stable (BIG) ★** |
| P35 | BL3 | EO25-A_20 | both years | growth |
| P36 | BL3 | EO25-B_21 | both years | **stable (SMALL) ★** |
| P37 | BL4 | EO38_1 | both years | ambiguous |
| P38 | BL4 | EO52_3 | both years | crash |
| P39 | BL4 | EO118_4 | both years | stable |
| P40 | BL4 | EO76_2 | both years | crash |
| P41 | BL4 | EO76_2 | 2025 only | 2025_only |
| P42 | BL5 | EO70_26 | both years | stable |
| P43 | BL5 | EO68-3_5 | both years | stable |
| P44 | BL5 | EO69_41 | 2026 only | 2026_only |

★ = recommended Next-test pair (BL3 within-BL stable
small-vs-big, D1 in the Executive summary and § Recommended field
pilot). SMALL / BIG tags in **bold** mark the three within-BL
stable pairs.

---

## Bottleneck Lineages (BL) — clustering and population naming (Step 30h)

**Why cluster populations into Bottleneck Lineages.** Phase 4
established a five-group Bottleneck Lineage (BL) partition at the
*EO* level (external `LEPA_EO_spatial_clustering` repo). The framework
reproduces that work at the *population* level, now with the two
years of pooled data instead of a single-year EO snapshot — a
check of how well the five-lineage grouping holds once populations
replace EOs.

**Approach.** Ward's D2 hierarchical clustering on the pairwise
haversine distances between the 44 population centroids (centroid
= size-weighted mean of event coordinates across 2025 + 2026).
Silhouette coefficient computed on the distance matrix for
*k* = 2 … 10.

**Result — five lineages, kept for continuity.** The silhouette
rises sharply from *k* = 4 (0.641) to *k* = 5 (0.716) and then
plateaus up to *k* = 9 (0.716 – 0.747), with its maximum at *k* = 8
(0.747). We keep ***k* = 5**, the number of lineages found by the
earlier EO-level analysis: the optimum is only 0.03 higher, and
*k* = 6 … 9 mostly subdivide the five lineages. *k* = 5 is fixed in
`step03_BL_clustering.py`, so BL identities stay stable as data are
added.

<a id="fig-3"></a>
![Figure 3](figures/step30h_dendrogram.png)

**Figure 3.** Ward's D2 dendrogram of the 44 populations
on pairwise haversine distances (y-axis = linkage distance in
metres). The dashed horizontal cut at *k* = 5 defines the five
Bottleneck Lineages used throughout the framework. Leaf labels are the
*post-numbering* populationIDs (populationID = position within
the BL on this dendrogram, numbered sequentially from BL1 to BL5)
so that adjacent IDs within a BL are also adjacent on the
clustering tree. Source: `step03_BL_clustering.py`.

<a id="fig-4"></a>
![Figure 4](figures/step30h_silhouette_curve.png)

**Figure 4.** Silhouette curve for the Ward's D2 clustering,
*k* = 2 … 10. The silhouette is highest at *k* = 8 (0.747), on a
plateau from *k* = 5 to 9 (0.716 – 0.747). We keep **five
lineages** (*k* = 5, silhouette 0.716), the number found by the
earlier EO-level analysis, for continuity; the optimum is only 0.03
higher. Source: `step03_BL_clustering.py`.

**BL definition rule — area DESC → connectivity DESC.** The five
Ward clusters are relabeled **BL1 … BL5** by (i) convex-hull area
of the member populations (DESC) and (ii) fraction of within-BL
pairs ≤ 5 km (DESC) — the same ordering convention used in the
external LEPA_EO_spatial_clustering repo so the BLs
remain comparable with the earlier EO-level labels. **Area is the
primary N<sub>e</sub> proxy** (feedback locked 2026-05-21).

| BL | Populations | Convex hull | Within-BL ≤ 5 km | Total `n_fertile` |
|----|:-:|:-:|:-:|:-:|
| **BL1** | 17 | 208.6 km² | 26 % | 4854 |
| **BL2** | 9 | 119.3 km² | 50 % | 5015 |
| **BL3** | 10 | 51.2 km² | 29 % | 4591 |
| **BL4** | 5 | 47.1 km² | 10 % | 1029 |
| **BL5** | 3 | 0.2 km² | 100 % | 592 |

BL1 is the largest-area, moderately-dispersed grouping; BL5 is
the smallest, most-compact grouping (every within-BL pair sits
within 5 km). The drift-collapsed EO24 group (P30 – P33) sits in BL3.

**Within-BL numbering: dendrogram leaf order.** Within each BL,
populations are renumbered **sequentially by their dendrogram
leaf position** (left → right on Figure 3). The populationIDs
BL1 holds (P1 … P17) are therefore mutually adjacent on the tree,
BL2 holds P18 … P26, and so on. This replaces the pre-BL
"populationID assigned in creation order" numbering and means
that adjacent populationIDs within a BL are also spatial-
clustering neighbours — making selection of comparable within-BL
pairs for the Objective 3 test design (§ C.5–C.6) direct from
the row ordering of every per-population figure.

<a id="fig-5"></a>
![Figure 5](figures/step30h_overview_by_BL.png)

**Figure 5.** The 44 populations across the Snake River
Plain, coloured by Bottleneck Lineage. Dot size ∝ total
`n_fertile` across 2025 + 2026; convex hulls drawn per BL
in matching colour. BL1 (17 populations, purple, N-E group,
largest area), BL2 (9 populations, blue, S group), BL3 (10
populations, red, central-W group, including the drift-collapsed
EO24 group), BL4 (5 populations, orange, NW group), BL5 (3
populations, green, far-W, most compact). Source: `step09_BL_overview_map.py`.

**Detailed maps per Bottleneck Lineage (local only).** The same step
draws one map per BL: element occurrences as convex hulls, populations
as polygons, and one zoom panel per population with its 50 m demes in
each season (2025 dashed, 2026 solid). They show slick-spot positions
in detail, so they are kept off the public repository: running
`./run_pipeline.sh` writes them to `figures/step30h_BL_map_BL1…5.png`
and `figures/step30h_BL_maps.pdf` on the local machine.

---

## Across-year population dynamic (Step 30g)

**Why classify populations by trend.** Once populations are
defined (§ Populations) and BL-grouped (§ above), the next
question is *how each population is faring across 2025 and 2026*.
That trajectory is what tells us which populations are strong
sampling candidates for the Objective 3 test design: large and
stable populations supply the "high-diversity" arm; crashing
populations supply the "drift-collapsed" arm; growing populations
diagnose which patches respond to good years.

**Trend classes.** Each population is tagged by its 2025 → 2026
`N_fertile` trajectory and diversity + pollen-compatibility
95 % CI overlap, under the equal-visitation-effort assumption
(both years received the same site-visitation protocol; the extra
2026 occurrence records were collected for seed banking and are
not used in this classification). The six classes:

| Class | Definition | n |
|---|---|:-:|
| **crash** | both-year; `N_fert_2026 ≤ 25 % of N_fert_2025`; Δdiversity < 0 | 3 |
| **stable** | both-year; 95 % CIs overlap on both metrics; `N_fert_total ≥ 20` | 15 |
| **growth** | both-year; `N_fert_2026 ≥ 4 × N_fert_2025` | 5 |
| **ambiguous** | both-year; none of the above | 3 |
| **2025_only** | present in 2025, absent above-ground in 2026 | 7 |
| **2026_only** | absent above-ground in 2025, present in 2026 | 11 |

**0 of 26 both-year populations show a pollen-compatibility CI
mismatch** between years — the structural redundancy of
§ C.0.b (Class I dominance + empirical zygosity) dominates. All
cross-year signal therefore lives on the **SRK-diversity axis**.

**Standout crash candidate.** `populationID 40 = locationCode EO76`
(one of the three § C.0.a clean-overlap EOs) collapsed from
**464 → 23 plants**. If the 2026 genotyping confirms the same
severely-collapsed allele pool observed at EO76 in 2025 (9 / 32
Fgs), that is **direct second-year evidence for H2** — the
diversity gap is persistent, not a one-year sampling fluke.

**Preliminary LARGE + SMALL candidate pair for the test design.**

- **LARGE = populationID 16** (EO27-1 + EO27RT merged), N_fert
  549 → 1442, predicted pollen compatibility 0.79 / 0.79,
  diversity 29 / 32.
- **SMALL = populationID 6** (EO67), N_fert 10 → 12, predicted
  pollen compatibility 0.71 / 0.71, diversity 7 / 7 (stable
  at a tiny scale).

<a id="fig-6"></a>
![Figure 6](figures/step30g_populations_classified.png)

**Figure 6.** The 44 populations grouped by across-year
trend class, each class as a vertical panel. Within each class,
rows are sorted by **BL (BL1 → BL5)** then **populationID**; the
BL grouping is shown as a left-margin **facet strip** (thin
coloured stripe + vertical BL label in BL colour) so BL blocks
read as sub-facets inside each class panel. Row labels are plain
short `P{N}` identifiers; the full `populationID → locationCode`
map is in
[`step30g_populations_classified.tsv`](tables/step30g_populations_classified.tsv).
Blue bars = `N_fertile_2025`, orange bars = `N_fertile_2026`.
Thin horizontal separators mark BL boundaries within each class.
**Headline patterns:** the three crash populations are in BL4
(P38, P40) and BL1 (P10); the fifteen stable populations span all
five BLs (BL1 5, BL2 4, BL3 3, BL5 2, BL4 1); the five growth
populations are in BL3 (3) and BL2 (2); ten of the eleven 2026-only
populations are BL1, mostly EO27 variants that emerged in 2026 after
a dry 2025. Source: `step06_populations_classified.py`.

<a id="fig-7"></a>
![Figure 7](figures/step30g_across_year_scatter.png)

**Figure 7.** Across-year prediction scatter per population
(both-year populations only, n = 26). **Panel A** — predicted
SRK diversity, 2025 (x) vs 2026 (y); **Panel B** — predicted
pollen compatibility, 2025 (x) vs 2026 (y). Dot size ∝ total
`n_fertile`; purple = both 95 % CIs overlap (stable); vermillion
= CI mismatch (prediction shifted). Equality diagonal shown.
Panel B's complete absence of CI mismatches (0 / 26) reaffirms
§ C.0.b: pollen compatibility is robust year-to-year. Panel A's
diversity CI mismatches concentrate the H2 signal (per-deme drift
beyond the species-wide prior) — see § C.0.c. Source:
`step11_year_scatter.py`.

---

## Delineating the deme — pollen-flight radius (Step 29b)

### The deme as a working hypothesis

Everything the pipeline predicts rests on how a **deme** is
delineated: **SRK allele diversity is a per-deme count** (how many
Fgs physically live in it), and **pollen compatibility is a
per-deme rate** (what fraction of pollen × stigma combinations
succeed inside it). Both aggregate to the location by set union and
size-weighted mean respectively (Figures 9, 10). **Mate limitation
at a location — the thing we actually want to test — is therefore
inherited directly from the deme-level predictions. Get the deme
wrong and the whole causal chain is wrong, so delineating demes
is paramount.**

The theoretical scaffolding is Wright's genetic-neighbourhood
concept (Wright 1943, 1946; Levin & Kerster 1974; Vekemans & Hardy
2004), which identifies the pollen-flight scale as the biologically
meaningful partitioning scale. We do **not** estimate Wright's
neighbourhood parameter *N*<sub>b</sub> itself — that would require
parent-offspring dispersal distances or a fine-scale *F*<sub>ST</sub> ~ distance
curve, neither of which exists for LEPA. Instead we delineate
**operational demes** using a hard-threshold connectivity rule
(50 m, bridged by events), calibrated by the sensitivity sweep
below. Our partition is a **geographic / topological upper bound**
on realised gene flow — flow across 50 m could be lower, but not
higher.

**The deme definition is a testable working hypothesis.** The
Part C seed-genotyping design compares the per-population pollen-compatibility prediction against observed per-mother compatibility.
A systematic location-level mismatch is not a framework failure —
it is information about **how good the 50 m deme definition
actually is**, and a future phase with parentage data or fine-scale
kinship markers could then refine the deme (narrower σ², extended
radius, behavioural weights) to close the gap. The deme is both
the hypothesis the pipeline *uses* to generate predictions and the
hypothesis the Part C data will *test*.

### Picking the pollen-flight radius

**Question.** With the deme defined operationally as "events whose
plants are connected within one pollinator flight", what *radius*
should the flight be? The answer must come from pollinator biology,
not from the sampling budget, so we sweep candidates and look for
the biologically defensible plateau.

**Approach.** Sweep candidate radii from 10 m to 200 m. At each
radius, build a within-location graph in which two events are
connected if any pair of their plants sits within that radius, and
record (i) the fraction of adults in a multi-event pollen-flow
component (connectivity), (ii) the fragmentation-aware sampling
cost ([`Sampling_design.md`](Sampling_design.md)), and (iii) the predicted random-mating
pollen compatibility under the sporophytic Class I/II model.

**Result.** **50 m is chosen on pollinator biology** (halictid /
small-bee foraging range), not on cost. The sweep shows how much the
design depends on that choice. At 50 m the design (mothers for 90 %
SRK allele detection per deme, [`Sampling_design.md`](Sampling_design.md)) needs **492 mothers
in 2025 (827 in 2026)**, 19 % (27 %) fewer than at 10 m. Assuming a
radius beyond 75 m would change it little (467 → 433 in 2025). If the
true radius is only 25 m, the 50 m design is **11 % short** (544 vs
492 mothers in 2025). Predicted pollen compatibility
does not depend on the radius (§ 30.2). Connectivity plateaus at
≥ 75 m, so 50 m still keeps meaningful fragmentation variation across
populations and Bottleneck Lineages.
Every downstream metric in this doc — SRK diversity prediction
(Figure 9), pollen compatibility prediction (Figure 11) — is
computed at this 50 m choice.

<a id="fig-1"></a>
![Figure 1](figures/step30_A_radius_sensitivity.png)

**Figure 1.** Pollinator-radius sensitivity sweep, **2025 and 2026 overlaid**. Four panels showing how connectivity, sampling cost (mothers needed for 90 % SRK allele detection per deme), predicted pollen compatibility, and the sustainable-band fraction of populations change across radii from 10 to 200 m. Each population is one unit; the radius changes only how its events split into components. 2025 draw = open circle + dashed line (33 populations); 2026 draw = filled circle + solid line (37 populations). Both years agree on the key qualitative shape: connectivity (Panel A, median population) rises steeply from 10 m, reaches 0.91 (2025) and 0.90 (2026) at 50 m and plateaus at 75 – 100 m; sampling cost (Panel B) drops steeply then flattens around 75–100 m; pollen compatibility (Panel C) is **radius-independent** at ≈ 0.77 under empirical zygosity in both years; the sustainable-band fraction (Panel D) sits at ≥ 97 % at every radius in both years. **50 m** is the primary radius adopted across the pipeline, chosen on pollinator biology (halictid / small-bee foraging range), not on cost: Panel B shows how sensitive the design is to that choice — 24 % (2025) / 32 % (2026) fewer mothers than at 10 m, little change beyond 75 m, and an 11 % shortfall if the true radius were only 25 m. Source: `step20_pollinator_radius.py`.

---

## Deme structure per population — what the 50 m partition yields (Step 29c)

### How the deme is built

**Question.** Within a population, how many demes does a plant
belong to, and how big is each one? A **50 m connected component**
is the **deme**: the set of adult plants whose events
are reachable from each other at the primary pollinator radius.
Plants in the same deme share pollen; plants in different
demes of the same population do not. The deme is the
**drift unit** on which every Phase A prediction is built.

**Approach.** For each population and year, build an event-level graph in which
two events are connected if any pair of their plants sits within
≤ 50 m. Connected components of that graph are the population's
demes. For each deme *c*, `component_N_fertile_c` is
the sum of `n_fertile_e` across its events — the adult count that
drives the per-component drift simulations in § 30.1 and § 30.2.
Reference: `step07_event_to_deme.py` (via `srk_units.py`) writes the
event → deme lookup (`step29a_event_to_deme.tsv`);
`step18_mating_pool_structure.py` builds the display here.

**Result (dataset-wide — how many demes do the 2025 populations hold?).**
The connectivity rule yields **99 operational demes across 33
populations** in 2025 at 50 m. The per-population breakdown:

| Demes in a population | Populations | Interpretation |
|:---:|:---:|:---|
| 1 | **13** | single connected deme; population-wide census = breeding unit |
| 2 | 9 | 2 sub-demes separated by > 50 m internal gaps |
| 3 | 3 | 3 sub-demes separated by > 50 m internal gaps |
| 4 | 3 | 4 sub-demes separated by > 50 m internal gaps |
| 5 | 1 | 5 sub-demes separated by > 50 m internal gaps |
| 6 | 1 | 6 sub-demes separated by > 50 m internal gaps |
| 9 | 1 | 9 sub-demes separated by > 50 m internal gaps |
| 11 | 1 | 11 sub-demes separated by > 50 m internal gaps |
| 16 | 1 | 16 sub-demes separated by > 50 m internal gaps |

- **13 / 33 populations operate as a single population-wide
  deme**, and **20 / 33 require within-population subdivision**
  into 2 – 16 operational demes to represent gene flow correctly.
- Deme sizes span **1 adult (SI floor) → 434 adults**; median
  **24 adults**.
- **29 of 99 demes (29 %) sit below the 8-plant species-pool
  threshold** (4 alleles per tetraploid plant × 8 plants = 32
  copies — the minimum needed to physically carry one copy of every
  SRK allele in the species pool). These demes are drift-limited
  for the 32-allele species pool regardless of location mean.
- **6 of 99 demes (6 %) sit at the N = 1 single-plant SI floor** —
  a single plant has nobody to mate with at the 50 m radius.
- Pattern across BLs: BL2 and BL3 carry the small-deme burden (11 of
  31 and 12 of 30 demes below 8 plants); BL3 and BL1 hold the most
  subdivided populations (P34 EO18-7/8: 16 demes; P16 EO27-1/RT: 11);
  BL4 is mostly fragmentation-robust (1 of 10 demes below 8); the two
  BL5 populations sampled in 2025 are one deme each.

**Why this matters.** Every Phase A prediction below is built on
this structure: the per-component SRK diversity prediction (Figure 9)
simulates each deme independently and unions the Fg sets at
the location level; the per-component pollen compatibility prediction
(Figure 11) runs the sporophytic simulation per deme and reports a
size-weighted mean per location; the fragmentation-aware sampling
allocation ([`Sampling_design.md`](Sampling_design.md)) assigns mothers per deme, with the
≥ 1-mother-per-event floor layered on top. The per-deme view is
the fastest way to anticipate which populations will stand out in
Figures 9 and 11.

**Per-population per-year deme-size view.** The dataset-wide
census above is a single-year snapshot. For a given population,
the number of demes **and their individual sizes** can shift
year-to-year because the above-ground set of occupied slick spots
flickers (annual plant, seed-bank-mediated emergence). Figure 8
drills down from the single-number per-population deme count
into the full per-year deme-size distribution — essential for
spotting populations whose predicted diversity shifts between
years not because the species-wide prior changed but because the
deme partition itself did.

<a id="fig-8"></a>
![Figure 8](figures/step29a_demes_per_population_year.png)

**Figure 8.** Per-population per-year deme-size distribution. One
row per population, panels stacked by BL (BL1 → BL5). Each row
holds up to two parallel strips of dots: **open circles = 2025
demes, filled circles = 2026 demes** (same open/filled convention
as Figures 9 and 11); one dot per deme at x = `component_N_fertile`
(log-scale deme size). Three vertical reference lines: red at
**N = 1** (single-plant SI floor — a lone plant has nobody to
mate with at 50 m), grey at **N = 8** (species-pool floor —
4 alleles per tetraploid × 8 plants = 32 copies, the physical
minimum to carry the full 32-Fg species pool), dashed at **N = 32**
(species ceiling of distinct Fgs). Left margin: short `P{N}`
identifier; right margin: per-year `{K demes}d / {N adults}`.
A row with dots in only one strip is a single-year population
(2025_only or 2026_only); rows where the 2025 and 2026 strips
sit at different x positions reveal populations whose deme
partition shifted year-to-year even where the population persists
across years. Source: `step10_demes_per_population_year.py`.
Data: [`step29a_demes_per_population_year.tsv`](tables/step29a_demes_per_population_year.tsv).

## Step 30 Phase A — Per-population predictions

### 30.1 Predicted SRK diversity per population

**Question.** How many SRK alleles does each population physically
hold under drift? Driven **per deme** by its fertile plants, then
unioned to the population, this is the drift output of the causal chain
and the raw material on which the pollen-compatibility prediction
(§ 30.2) operates. Whether our sampling can recover these alleles is
covered in [`Sampling_design.md`](Sampling_design.md) § 3.

**Approach — how the simulation works, step by step.**

Think of the species prior P1 as a bag with 32 kinds of coloured
balls — one colour per Fg. The bag is **not evenly filled**: FG001
occupies 41 % of it, the next five Fgs together another ~25 %, and
the 20+ rare Fgs each fill under 2 %. Those shares are the
*species-wide empirical frequencies* measured from the Canu-amplicon
263-individual inventory. Each LEPA location only gets a limited
handful of draws from that bag, so small locations are almost
guaranteed to miss the rare colours by chance — that is **genetic
drift**. The prediction turns this intuition into numbers in three
steps, run for every 50 m deme:

1. **Build the component's local pool (drift step).** For a
   component with `component_N_fertile_c` fertile plants, draw
   `4 × component_N_fertile_c` balls from P1 (because every plant is
   tetraploid and carries 4 SRK alleles). The unique colours drawn
   are the Fgs **present at that component**; the proportions in the
   draw are the component's local Fg frequency vector `f_local_c`.
   Small `N_c` → few draws → missing some rarer Fgs.

2. **Sum components to the location (union step).** The location
   holds a Fg if *any* of its components holds it, so the location's
   pool size on this replicate = |union of per-component present
   sets|. That is why diversity is unioned (set operation), not
   averaged (continuous-rate operation).

3. **Repeat 1 000 times, report mean and 95 % credible interval.**
   Each replicate uses an independent draw from P1; the 2.5 % and
   97.5 % percentiles across replicates give the credible interval
   seen as error bars in Figure 9.

<a id="fig-9"></a>
![Figure 9](figures/step30h_pred_diversity.png)

**Figure 9.** Predicted SRK allele diversity per population, 2025 and 2026 overlaid. One row per population; the 2025 draw is shown as an open circle, the 2026 draw as a filled circle in the same BL colour, with a thin connector line joining the two when both years are present. 95 % credible intervals shown as horizontal error bars; dot size ∝ `N_fertile` for that year. Single-year populations show only the applicable year's marker and no connector. Dotted vertical line = species-wide ceiling (32 Fgs). Panels stacked by BL (BL1 → BL5, area DESC → connectivity DESC); within each BL populations ordered by the mean of their available-year predicted means. Left margin: short `P{N}` identifiers; right margin: per-year `{K demes}d / {N adults}`. The full `populationID → locationCode` map is in [`step30g_populations_classified.tsv`](tables/step30g_populations_classified.tsv). The across-year comparison within a population is the within-pipeline H2 test of § C.0.c — populations whose predicted diversity shifts markedly between years carry the strongest signal of per-deme drift beyond the species-wide prior. Source: `step12_predictions_by_BL_year.py`.

### 30.2 Predicted pollen compatibility per population

**Question.** Under the sporophytic Class I / II + empirical-zygosity
model, what fraction of pollen would a mother at each location be
compatible with under random mating?

**The underlying biology — sporophytic Class I / II dominance.**
LEPA SRK alleles fall into two dominance classes. **Class I is
dominant**: a plant carrying any Class I allele expresses only its
Class I alleles on both pollen and stigma, silencing its Class II
alleles. A plant carrying only Class II alleles co-expresses all
of them. **A cross is rejected if the two parents share any
expressed allele**; between-class crosses (Class I plant × Class II
plant) are always compatible by construction. This single rule
drives the whole § 30.2 prediction — and crucially, it is what
makes the § C.0.b hypothesis decomposition work: without Class I
dominance, within-class drift concentration (EO70's FG024) and
zygosity composition (EO76's homozygous buffer) would be
indistinguishable from a plain allele-frequency effect. Figure 9
walks through the rule in three panels before we get into the
simulation steps.

<a id="fig-10"></a>
![Figure 10](figures/step30_A_si_model_schematic.png)

**Figure 10.** Sporophytic SI with Class I / Class II dominance in tetraploid LEPA. **Panel A — dominance within one plant.** Case A: a plant with ≥ 1 Class I allele expresses only its Class I alleles; its Class II alleles are silent (shown faded). Case B: a plant carrying only Class II alleles expresses all four Class II alleles co-dominantly. **Panel B — between-plant recognition, worked example.** Mother M carries {FG001, FG002, FG024, FG031}; her Case-A expressed set is {FG001, FG002}. Three candidate fathers: F1 shares FG001 with M → rejected; F2 is all-Class-II so between-class → always compatible; F3 shares FG002 with M → rejected. **Panel C — compatibility rule by cross type.** Class I × Class I: compatible if their expressed Class I alleles differ (shared Class II is irrelevant because Class II is silent on both sides). Class I × Class II: always compatible by construction (disjoint expressed classes). Class II × Class II: all four alleles expressed on both sides, compatible only if none are shared.

**Approach — how the simulation works, step by step.**

Same per-component logic as § 30.1, but the recognition rule changes
from "did we observe this Fg at least once?" to "if this mother's
expressed Fgs are {a, b}, what fraction of pollen in her component
would she accept?" Five steps per replicate:

1. **Build the component's local pool (same as § 30.1 step 1).**
   Draw `4 × component_N_fertile_c` alleles from P1 → local Fg
   frequencies `f_local_c`. This is the drift-collapsed pollen pool
   the component's mothers see.

2. **Sample mother genotypes under the empirical LEPA zygosity.**
   Each mother gets 4 SRK alleles drawn from `f_local_c`, but
   combined into tetraploid genotypes that match the observed LEPA
   mix: 66 % of mothers end up with 1 distinct Fg (homozygotes),
   32 % with 2 distinct, 2 % with 3 distinct (from
   `srk_zygosity_empirical.tsv` — see Part C § C.0 for how this was
   measured from the 367 Canu-amplicon adults).

3. **Compute per-mother pollen compatibility under sporophytic Class
   I / II SI.** A pollen parent expresses only Class I alleles if it
   carries any, else all its Class II alleles co-dominantly. A cross
   is rejected if parents share any expressed allele. For a mother
   expressing Fg set *M*:
   - Mother has ≥ 1 Class I allele: `P_compat = (1 − p(M))⁴` where
     p(M) is the local frequency sum of her expressed Fgs.
   - Mother is pure Class II: `P_compat = 1 − (1 − p_I)⁴ +
     (1 − p_I − p(M))⁴`.

   Average across the mothers sampled at this component → component
   P_compat on this replicate.

4. **Aggregate components to the location (weighted mean, not
   union).** P_compat is a continuous rate, so a location's number
   is the component sizes weighted mean:
   P_compat_location = Σ_c (P_compat_c · component_N_fertile_c) /
   Σ_c component_N_fertile_c. Big components dominate; small
   components still count — that is where drift pulls the location
   mean down.

5. **Repeat 400 times, report mean and 95 % credible interval.**

Why no sampling term enters P_compat: § 30.2 asks "what fraction of
pollen is compatible with a random mother at this component?" — a
biological property of the pollen pool, not a function of how many
seeds we sequence. Seed counts enter § 30.1 (sampling detection) and
Phase B (the regression against observed seed set), not here.

Both tables are emitted:
[`step30_A_prediction_location_pcompat.tsv`](tables/step30_A_prediction_location_pcompat.tsv)
(headline per location) +
[`step30_A_prediction_component_pcompat.tsv`](tables/step30_A_prediction_component_pcompat.tsv)
(one row per component, so struggling sub-components inside a
sustainable-mean location remain visible).

**Result.**

- **Species-mean P_compat = 0.78** (sustainable band by construction).
- **Traffic-light bands** (1/3 and 2/3 of species mean): failed
  < 0.23, struggling 0.26–0.52, sustainable ≥ 0.52.
- **The EO24 group (P30 – P33, BL3)** is the main conservation
  concern — predicted pollen compatibility 0.13 – 0.51 at locations of
  1 – 3 plants, the lowest in the dataset (EO24-2, one plant: 0.13;
  EO24 and EO24-1: 0.37 – 0.39; EO24-7: 0.51), with credible intervals
  reaching down to 0. The smallest EO26-3 location (5 plants) follows
  at 0.49.
- Every other location sits in the sustainable band at the mean; the
  tightest credible intervals belong to the largest locations
  (EO27-1, EO32).

**Connection to Figure 9.** Uses the same drift mechanism as the
diversity figure, except drift is now modelled **per component**
(one pool per 50 m connected component sized by
`4 × component_N_fertile`) rather than through a single
"largest-component" proxy. Figure 9 shows that the existing LEPA
seed data recovers the true local Fg pool at every location — that
validation transfers directly to Figure 11, meaning the location-
mean P_compat from the per-component simulation is a faithful
estimator of the population-mean P_compat.

<a id="fig-11"></a>
![Figure 11](figures/step30h_pred_pcompat.png)

**Figure 11.** Predicted per-population pollen compatibility under sporophytic Class I / II SI + empirical LEPA zygosity, 2025 and 2026 overlaid. Same layout as Figure 9: open circle = 2025, filled circle = 2026, same BL colour, thin connector line for both-year populations. 95 % CI horizontal error bars; dot size ∝ `N_fertile` for that year. Traffic-light bands: red = failed (< 0.259), amber = struggling (0.259 – 0.519), green = sustainable (≥ 0.519). Dotted green line = sporophytic species mean 0.778. Panels stacked by BL (BL1 → BL5); within each BL rows ordered by the mean of available-year means. Left margin: short `P{N}` identifiers; right margin: per-year `{K demes}d / {N adults}`. The full `populationID → locationCode` map is in [`step30g_populations_classified.tsv`](tables/step30g_populations_classified.tsv). Nearly all populations with ≥ 8 adults sit on or just below the species mean — the structural redundancy of Class I dominance + empirical zygosity (§ C.0.b) swamps between-population variation; BL3's EO24 tail (P30 – P33, 1 – 3 plants each) is the only group with meaningful credible-interval spread below the mean. All 26 both-year populations have overlapping 2025 and 2026 CIs on pollen compatibility (0 / 26 mismatches), so the across-year differences are concentrated in the single-year populations. Source: `step12_predictions_by_BL_year.py`.


### 30.3 Do populations produce more or fewer seeds than their size predicts? Mate limitation vs SI escape across all 44 populations (Step 30i)

**Why this test exists.** The SRK-based Part C validation tests (§ C.0.a – § C.0.c, below) are confined to the six genotyped sites (EO67 = P6, EO70 = P42, EO76 = P40 + P41, and P34, P13, P16 of EO18 / EO27) because they need SRK genotypes. **This step adds a behavioural check that uses no SRK data** and therefore runs on **all 44 populations NOW**, before any seed-genotyping campaign starts — a direct validation of the Phase A pollen-compatibility prediction that places § 30.3 cleanly inside the predictions pipeline rather than waiting for Part C. It rests on one hypothesis: **seed production reflects the reproductive status of the population** — under-production signals mate limitation, over-production signals SI escape.

**Approach — across-species build, then within-population test.**
The test splits into two stages:

1. **Across-species build**: Stage 1 picks the best size
   predictor out-of-sample via 10-fold CV RMSE on log<sub>10</sub>
   seed yield across four candidate allometries — `height`,
   `crown`, `area = π/4 · crown · height`, and `crown + height`
   (free exponents). Stage 2 refits the CV-winner on **all plants
   pooled across populations and years**: this is the
   species-wide expectation curve. Every plant then has an
   expected log<sub>10</sub>(yield) and a residual = observed −
   expected. **Figure 12a** shows this build: observed vs
   expected for every plant, log-log, coloured by BL, with the
   1:1 species curve and the ×½ / ×2 reference lines.
2. **Within-population test**: Stage 3 tests per (populationID,
   year) whether the mean residual departs from 0 (< 0 → mate
   limitation, > 0 → SI escape) using
   a one-sample t-test (n ≥ 3) with a Benjamini-Hochberg FDR
   across the ~ 80 strata. **No min-n threshold** — small
   populations are **included and highlighted** (user focus
   2026-10-07): small-N rows carry wide CIs, but still show up
   in the plot. **Figure 12b** is the per-population forest
   plot of fold-of-expectation (= 10<sup>mean residual</sup>).

Each plant is its own control via the species-wide allometry, so
a population that shows a systematic residual is NOT just a size
artefact — it is behaving differently from same-size plants
elsewhere in the species.

**Caveat on 2026 germplasm coverage.** A fraction of the 2026
mother-plant occurrences do not yet have associated germplasm
records attached in the DB (clean-up in progress 2026-10-07).
2026 shortfall calls may therefore be inflated at populations
that are still missing collections; re-run step30i once the
2026 germplasm cleaning is complete.

<a id="fig-12a"></a>
![Figure 12a](figures/step30i_species_calibration.png)

**Figure 12a — the across-species build.** Species-wide
calibration scatter underlying the per-population test. **Panel
A (left):** every plant that has crown + height + seed yield
plotted as expected (species-wide allometry, log scale) vs
observed (log scale), coloured by its population's BL. Dashed
line = the species curve (fold = 1); dotted lines = ×½ and ×2
reference bands. **Panel B (right):** the model build itself —
the Stage 1 10-fold CV predictor-selection table (which
allometry wins out-of-sample) + the Stage 2 refitted
expectation-model summary (coefficients, R², residual SD). This
is what every row of Figure 12b is measured against. Source:
`step13_size_seed_yield.py`.
Data: [`step30i_size_seed_cv_rmse.tsv`](tables/step30i_size_seed_cv_rmse.tsv).

<a id="fig-12b"></a>
![Figure 12b](figures/step30i_size_seed_population.png)

**Figure 12b — the within-population test.** Observed seed yield
per plant divided by the expected yield from Figure 12a's
species-wide allometry, aggregated per (populationID, year). One
row per population; **open circle = 2025 draw, filled circle =
2026 draw** in the same BL colour, slightly offset vertically;
dot size ∝ √*n* plants. Horizontal bars = 95 % credible
interval; the dashed vertical line at fold = 1 is the species
curve. The light grey ×½ and ×2 reference bands are the
operational "biologically meaningful shortfall / surplus"
markers. **Two halo colours flag the FDR-significant
predictions**, both tied to the specific (populationID,
**year**) data point and both informative tests of the
pollen-compatibility prediction:

- **Red halo = below expectation** (q < 0.05 and fold < 1) —
  the population's plants produced fewer seeds than their sizes
  predict, a behavioural signature of **mate limitation
  surfacing**.
- **Green halo = above expectation** (q < 0.05 and fold > 1) —
  the population's plants produced more seeds than their sizes
  predict, a behavioural signature of **SI escape**:
  self-compatible plants setting self-seed on top of what
  outcrossing supplies.

The flag is per-year, not per-population: a population halo'd
in only 2025 means the 2025 draw was flagged and 2026 was not
(or vice versa); a population halo'd in both years means the
same signal is sustained across years. Panels stacked by BL;
left margin = `P{N}`; right margin = per-year `n` plants.
Source: `step13_size_seed_yield.py`. Data:
[`step30i_size_seed_population_strata.tsv`](tables/step30i_size_seed_population_strata.tsv).

**Framing: mate limitation vs SI escape.** The test separates
the 44 populations into **reproductive-status buckets** by
whether their observed seed yield falls below, at, or above
what plant size alone predicts:

- **Mate-limited populations** (red halos in Figure 12b) —
  plants produce fewer seeds than their sizes predict. The
  pollen-compatibility pathway is surfacing behaviourally — the
  § C.0.b buffering stack is NOT absorbing the drift pressure
  at this population. **These are the conservation-priority
  populations to intervene on.**
- **Reproduction intact** (CI overlaps fold = 1) — plants
  produce the seed their size predicts. The pollen-compatibility
  prediction holds: the § C.0.b buffering mechanisms (Class I
  dominance + class imbalance + tetraploid homozygosity) are
  keeping outcrossing working. **No priority intervention
  needed.**
- **SI-escape candidates** (green halos in Figure 12b) — plants
  produce *more* seed than their size predicts. Over-production
  is the predicted signature of self-fertilisation by
  self-compatible (SC / pSI) plants. **Priority for seed-parent
  genotyping** to confirm selfing — SI escape erodes the SI
  system the conservation programme relies on.

The one-shot "which populations to prioritize" table is
[`step30i_conservation_priority.tsv`](tables/step30i_conservation_priority.tsv),
one row per populationID. The test is run per (population, year);
the two years are then combined into **one population-level
status**, ranked for restoration (`priority_rank`, 1 = act
first; `restoration_action` gives the action):

| Rank | Status | n (pooled curve / robust) | Rule | Action |
|:-:|---|:-:|---|---|
| 1 | **mate_limited** | **8 / 6** | FDR-sig under-producing in ≥ 1 year, never over | Intervene (cross-pollination / founder augmentation) |
| 2 | **switching** | **1 / 1** | FDR-sig under in one year AND over in the other | Intervene; genotype seed + mothers from both years |
| 3 | **candidate_mate_limited** | 3 / 3 | fold < 0.8 in ≥ 1 year, not FDR-sig | Watch; sample more mothers |
| 4 | **si_escape_candidate** | 12 / 7 | FDR-sig over-producing in ≥ 1 year, never under | Genotype seed + mothers to confirm selfing |
| 5 | **reproduction_intact** | 10 / 8 | every testable year at expectation | No action |
| — | **year_sensitive** | — / 9 | status differs between the pooled curve and the year-offset curve | Re-test with 2027 data; genotype |
| 6 | **insufficient_data** | 3 / 3 | no year with n ≥ 3 plants | Sample more plants |

**The size → seed curve, and the year effect.** The curve is fitted
on both years pooled (1 441 plants). A year effect estimated *within*
populations (population random intercept) is now significant: 2026
plants made 15 % fewer seeds for their size than 2025 plants in the
same populations (fold 0.85, p = 0.005;
[`step30i_year_effect_check.tsv`](tables/step30i_year_effect_check.tsv);
it was 0.90, p = 0.07 before the 2026-10-20 germplasm batch). With two
years we cannot tell whether this is a year effect (weather) or a real
population-wide decline, which is itself a mate-limitation signal. The
test is therefore run twice — with the pooled curve and with the
within-population year effect as a fixed offset on 2026 (curve
refitted) — and **only statuses that agree are ranked**; the rest are
**year-sensitive** (`robust_priority` and `priority_year_adjusted` in
[`step30i_conservation_priority.tsv`](tables/step30i_conservation_priority.tsv)).
A naive year term (without the population structure) is not used: it
would read the 2025-only BL2 over-producers as a "good 2025" (fold
0.65) and remove real population differences.

**The 6 robust mate-limited populations** (rank 1):

- **P10 EO27-3 (BL1)** — flagged **both 2025 and 2026** (fold
  0.44 / 0.02) — most severe; the EO27 size-restriction field
  note confirmed behaviourally.
- **P42 EO70 (BL5)** — flagged **both years** (fold 0.58 /
  0.56) — behavioural reproduction of the predicted
  pollen-compatibility shortfall at the one clean-overlap
  population where the SRK test flagged a gap. Independent
  corroboration without any SRK data.
- **P8 EO27 (BL1)** — flagged 2025 (fold 0.44).
- **P4 EO30-3 (BL1)** and **P1 EO112 (BL1)** — flagged 2026
  (fold 0.23 and 0.34; 16 and 7 plants). Both are late-season
  populations seen only in 2026.
- **P32 EO24 (BL3)** — flagged 2026 (fold 0.22, 5 plants).

**The switching population** (rank 2) — over-producing in 2025,
under-producing in 2026, under both curves:

- **P19 EO8_27+28 (BL2)** — fold 1.74 (2025, n = 37) → 0.60
  (2026, n = 137); 1.61 → 0.65 with the year offset. The swing is
  far larger than the within-population year effect (fold 0.85).

Over-production followed by under-production could reflect SI
escape followed by inbred progeny, but the seed bank argues against a
one-year link: only ≈ 6 % of a seed cohort germinates the next year
(Meyer, Quinney & Weaver 2006), so the 2026 plants come mostly from
older cohorts. Turnover of the emerging plants, or weather, are the
likelier causes of a swing this size. Seed-parent genotyping of
seed lots from both years is the direct test.

**The 7 robust SI-escape candidates** (rank 4 — seed set
significantly above expectation under both curves):

- **BL2 2025 is surplus-dominant**: **P18 EO8_29+47** (fold
  2.26), **P23 EO26-3_31+32** (2.19), **P24 EO26-1_30** (1.65) and
  **P25 EO61_38** (1.91) — plus P19 in 2025 (switching) and P26 EO29
  (1.49, year-sensitive). BL2 is the mid-area / high-connectivity
  grouping — a lineage-wide over-production signal, the strongest
  SI-escape candidate set and the first place to genotype seed for
  selfing.
- **BL4**: **P40 EO76_2** 2025 (2.22) and **P37 EO38_1** (1.18 →
  2.47 in 2026). P40 at EO76 is particularly striking: surviving
  2025 adults over-produced *before* the population collapsed to 23
  plants in 2026 — consistent with selfing in a pool already short
  of compatible mates (EO76 holds 9/32 Fgs, § C.0.a).
- **BL1**: **P2 EO30-1_13** (0.92 → 1.50 in 2026).

**The 9 year-sensitive populations** — status depends on whether the
2026 drop is treated as a year effect:

- Mate-limited with the pooled curve: **P16 EO27-1 + EO27RT** (BL1;
  0.39 in 2025, 1.10 in 2026 — *switching* with the offset) and
  **P36 EO25-B** (BL3; 0.52 in 2026 — *candidate* with the offset).
  P36 and P16 are the small side of D1 and the big side of D3.
- SI-escape with the pooled curve, intact with the offset (all from
  their 2025 surplus): **P26 EO29** (1.49), **P29 EO32** (1.56),
  **P34 EO18-7 + EO18-8** (1.28; the big side of D1), **P35 EO25-A**
  (1.44) and **P43 EO68-3** (1.32).
- Intact with the pooled curve, candidate mate-limited with the
  offset: **P13 EO27-1** and **P27 EO48** (both 0.83, 2025).

**Candidate mate-limited (sub-threshold)** — fold < 0.8 but not
FDR-sig, worth watching: **P6 EO67_39** (fold 0.54 / 0.54,
SMALL candidate from § C.0.a — needs more mothers to confirm
statistically), **P3 EO30-4_52** (0.55, 2026) and **P20
EO26-2_35** (0.70, 2025).

**Re-run once the 2026 germplasm records are complete** —
2026 calls (including the switching population) may shift.

**Over-production also appears at tiny N** — small populations
with seed yield above the species curve despite N ≤ 20:
P41 EO76_2 2025 n = 4 fold 2.43; P30 EO24-7_25 2025 n = 3
fold 2.51. Neither is FDR-significant, but small, allele-poor pools are exactly where
SI escape is expected to pay off. The seed-set data separate the
"small → mate-limited" cases from the "small → SI escape"
candidates.

### Does seed yield follow the genetic predictions?

**Question.** Do populations predicted to have fewer compatible mates
(§ 30.2) or fewer SRK alleles (§ 30.1) produce fewer seeds for their
size?

**Approach.** For every population-year with ≥ 3 plants (45), rank
correlation (Spearman) between the fold of expectation (Figure 12b)
and the predicted pollen compatibility, predicted SRK allele
diversity and fertile-plant census.

**Result.** Seed yield does not track predicted pollen compatibility
(ρ = 0.17, p = 0.27, 45 population-years). The only hint is in 2026
(ρ = 0.45, p = 0.08, 16 population-years); in 2025 there is none
(ρ = 0.09, p = 0.64). Predicted SRK diversity (ρ = 0.17, p = 0.27)
and census (ρ = 0.17, p = 0.28) do not correlate either
([`step30i_seed_yield_vs_predictions_tests.tsv`](tables/step30i_seed_yield_vs_predictions_tests.tsv)).
An earlier version of the drift model, which drew each plant's four
SRK copies independently, gave ρ = 0.38 (p = 0.01); simulating plants
with the observed zygosity (2026-10-23) removed it. Predicted
compatibility spans 0.52 – 0.79 across the tested population-years,
and **the mate-limited populations sit at the species mean** (P8,
P10, P42: 0.77 – 0.79). The drift model does not predict where mate limitation
occurs. The genotyped sites (§ C.0.a) suggest why: local demes have
lost far more SRK alleles than one draw from the species-wide
frequencies predicts. Seed yield is therefore the direct read-out of
reproductive status, and the mismatches are the first genotyping
targets.

<a id="fig-12c"></a>
![Figure 12c](figures/step30i_seed_yield_vs_predictions.png)

**Figure 12c.** Seed yield against the genetic predictions, one point
per population-year (open = 2025, filled = 2026; colour =
Bottleneck Lineage). **Panel A** — predicted pollen compatibility;
**Panel B** — predicted number of SRK alleles. The y-axis is observed
seed yield relative to the plant-size expectation (log scale; dotted
line = at expectation). Mate-limited populations (below the line)
sit at the same predicted compatibility as the over-producers.
Source: `step13_size_seed_yield.py`. Data:
[`step30i_seed_yield_vs_predictions.tsv`](tables/step30i_seed_yield_vs_predictions.tsv).

---

## Sampling design

How many seeds per mother, how many mothers per deme, which demes to
sample, and whether the material we already hold can test the
predictions are covered in [`Sampling_design.md`](Sampling_design.md): seedlings per mother and statistical
power (§ 1), mothers per deme (§ 2, Step 1), the sampling check per
population (§ 3), the selection of spatially separated demes for the
model populations (§ 4, Step 2) and the locked recipe with the lab
files (§ 5).

---

## Step 30 Part C § C.0 — Empirical validation of the pollen-compatibility model

**The genotyped plants.** The SRK genotypes come from plants grown ex
situ from wild-collected seed, plus some leaf tissue sampled from wild
plants (2023). In the LEPA database they are placeholder occurrences
(no event, location or date yet), so each plant is placed in a
population through the genotyping sample records
(`srk_genotype_units.py`): its seed lot → wild mother → her 2025 slick
spot, else the sample's location, else the tissue-bank position
(nearest slick spot within 500 m).

The validation runs on two kinds of site:

- **Three EO sites** — EO67 (= P6), EO70 (= P42) and EO76 (= P40 +
  P41), compared at location scale (§ C.0.a, § C.0.b, § C.0.c).
- **Three populations of the split EOs** (added 2026-10-21) — every
  population of EO18 and EO27 with ≥ 10 placed plants: **P34 (EO18)**,
  **P13 and P16 (EO27)**, compared with their 2025 population-level
  prediction (§ C.0.a). P8 (5 plants) and P10 (4) of EO27 have too few.
  **EO25** stays out: its plants carry no seed-lot or location link, so
  they cannot be split between P35 (EO25-A) and P36 (EO25-B).

Seedlings of the same wild mother are counted as independent plants
(no family correction), as for the three EO sites.

**Background — the split EOs as whole EOs** (from
`step24_EO_scale_validation.py`, observed fathers from local allele
frequencies vs predicted from the species-wide prior, same observed
mothers for both):

| EO | plants | Observed pollen compatibility (95 % CI) | Predicted (95 % CI) |
|---|---:|---|---|
| EO27 | 61 | 0.77 (0.73–0.80) | 0.80 (0.77–0.84) |
| EO25 | 51 | 0.77 (0.74–0.80) | 0.80 (0.77–0.84) |
| EO18 | 39 | 0.69 (0.64–0.74) | 0.75 (0.71–0.79) |

Pooled over the whole EO, EO27 and EO25 pass and EO18 shows a mild
shortfall. Split into populations (§ C.0.a), EO27 is not homogeneous:
P13 falls far below the prediction and P16 sits above it.

**Caveat.** The species-wide prior was built from these same plants,
so species-mean alignment is guaranteed; what this tests is
**robustness to local drift**, not absolute calibration.

### C.0.a — Part C anchor: observed vs predicted at the genotyped sites

**Question.** Does the prediction reproduce the SRK diversity and
pollen compatibility observed in the genotyped plants of each site?

**The sites.** EO sites: EO67 (P6, 37 plants), EO70 (P42, 74), EO76
(P40 + P41, 76). Populations of the split EOs: P34 (EO18, 32 plants),
P13 (EO27, 12), P16 (EO27, 15). In total **246 genotyped plants**.

**Approach.** For each site, run `step25_partC_clean_overlap.py`:
1. Observed per-plant pollen compatibility under sporophytic Class I /
   II dominance + empirical zygosity, with fathers drawn from the
   **observed local allele frequencies** (not the species prior).
2. Observed number of distinct SRK alleles in the genotyped plants.
3. **No-drift upper bound** at the same sample size:
   `E[distinct alleles | 4·n plants drawn from the species prior]`.
   A large gap below it isolates **drift**, not sampling.
4. Compare with the prediction — location-level for the EO sites
   (`step30_A_prediction_location_pcompat.tsv`,
   `_diversity.tsv`), population-level for 2025 for the split-EO
   populations (`step30g_prediction_population_year.tsv`).

**Result.**

| Site | Plants | Observed pollen compatibility (95 % CI) | Predicted (95 % CI) | Observed alleles | No-drift upper bound | Predicted alleles (95 % CI) |
|---|---:|---|---|---:|---:|---|
| **EO67 (P6)** | 37 | 0.697 (0.660–0.737) | 0.62 (0.39–0.74) | **7** of 32 | 19.2 | 6.5 (4.0–10.0) |
| **EO70 (P42)** | 74 | **0.532** (0.510–0.559) | **0.78** (0.74–0.83) | **6** of 32 | 23.8 | 22.5 (18.0–27.0) |
| **EO76 (P40 + P41)** | 76 | 0.722 (0.694–0.753) | 0.78 (0.75–0.81) | **9** of 32 | 24.0 | 29.4 (27.0–32.0) |
| **P34 (EO18)** | 32 | **0.671** (0.613–0.722) | **0.79** (0.75–0.83) | **7** of 32 | 18.2 | 29.2 (26.0–31.3) |
| **P13 (EO27)** | 12 | **0.436** (0.310–0.577) | **0.79** (0.74–0.84) | **5** of 32 | 11.8 | 27.8 (24.0–31.0) |
| **P16 (EO27)** | 15 | 0.843 (0.816–0.871) | 0.79 (0.74–0.83) | **8** of 32 | 13.2 | 29.4 (26.0–32.0) |

- **EO67 (P6)** — observed and predicted intervals overlap on pollen
  compatibility, and the observed allele count sits inside the
  predicted interval: **the model passes** at this small BL1 site.
- **EO70 (P42)** — observed compatibility 0.25 below the prediction
  (intervals do not overlap); 6 of 32 alleles against 23 predicted and
  a no-drift bound of 24: **massive drift collapse in a large
  population**.
- **EO76 (P40 + P41)** — compatibility marginally below the
  prediction; 9 of 32 alleles against 29 predicted: **severe diversity
  collapse**, yet compatibility holds because FG001 (Class II, 41 %
  species-wide) makes most plants compatible even with few alleles.
- **P34 (EO18)** — compatibility below the prediction (intervals just
  separate) with 7 of 32 alleles: the same pattern as EO70, milder.
  P34 is the big side of the D1 model pair.
- **P13 (EO27)** — the largest shortfall of all (0.44 vs 0.79), with
  only 5 alleles and FG001 (0.56) and FG002 (0.27) dominating: only
  4 % of its SRK copies are Class I, so almost every plant expresses
  Class II (a class-balance collapse, § C.0.b). In the seed-yield test
  P13 is intact with the pooled curve and a candidate with the year
  offset (year-sensitive, § 30.3). Few plants (12, from 7 mothers):
  treat as indicative.
- **P16 (EO27)** — compatibility **above** the prediction here (0.84
  vs 0.79): its 8 alleles are evenly spread, with rare alleles (FG005,
  FG017) at 0.10 each. This comparison draws fathers with the
  species-wide zygosity; with P16's own zygosity (53 % of plants
  express ≥ 3 identities) compatibility falls to 0.77, below the
  species-wide value (§ C.0.b). Few plants (15, from 7 mothers):
  indicative.

At five of the six sites the allele count is far below the prediction
(5 – 9 observed against 22 – 29 predicted; P6 matches), while compatibility follows
the prediction at two sites (P6, EO76), falls below it at three
(EO70, P34, P13) and exceeds it at one (P16) — the decomposition in § C.0.b
explains why the allele count alone does not set compatibility.

**Outputs.**
- [`step30_B_partC_clean_overlap_per_location.tsv`](tables/step30_B_partC_clean_overlap_per_location.tsv) — one row per site (`unit` = EO or population; `prediction_source` names the prediction used).
- [`step30_B_partC_clean_overlap_fg_frequencies.tsv`](tables/step30_B_partC_clean_overlap_fg_frequencies.tsv) — long form, one row per (site, allele) with observed and species-wide frequency.
- [`step30_B_partC_clean_overlap.png`](figures/step30_B_partC_clean_overlap.png) — three-panel figure.

<a id="fig-13"></a>
![Figure 13](figures/step30_B_partC_clean_overlap.png)

**Figure 13.** Part C anchor at the six SRK-genotyped sites — three EO sites (**P6 = EO67**, **P42 = EO70**, **P40 + P41 = EO76**, the last compared as one site) and three populations of the split EOs (**P34 of EO18**, **P13 and P16 of EO27**). Panel order follows the causal chain: diversity → pollen compatibility → per-allele drift fingerprint. **Panel A — SRK diversity.** Predicted (x) vs observed in the genotyped plants (y), one square per site with 95 % interval error bars; 1:1 diagonal and ceiling at 32 alleles. Only P6 sits near the diagonal; every other site sits far below it. **Panel B — pollen compatibility.** Predicted vs observed, with the traffic-light background (red = failed, amber = struggling, green = sustainable) and 1:1 diagonal. P6, P40 + P41 and P34 lie close to the diagonal or just below; P42 and P13 are clear outliers below, P16 above. **Panel C — per-allele drift residual** (observed − species-wide frequency) per site, alleles sorted by species-wide frequency (most common → rarest). **Green = enriched**, **red = depleted**, **× = absent** at the site. P42 shows FG001 (+21 %) and FG024 (+18 %) with 26 of 32 alleles absent; P13 is dominated by FG001 and FG002; P16 carries several normally rare alleles at raised frequency; P6's small-population signature raises the rare FG018 and FG023 — a founder effect rather than classical drift. Source: `step25_partC_clean_overlap.py`.

### C.0.b Diversity collapse → pollen compatibility: hypothesis decomposition

**Why this analysis exists.** Panel A of § C.0.a (SRK diversity
predicted vs observed) exposed a disconnect that pollen
compatibility alone would have hidden: EO70 (P42) and EO76 (P40 +
P41) lose 17 and 20 of their predicted alleles, but observed pollen
compatibility drops by 0.24 at P42 and barely 0.03 at EO76. P6, with
7 alleles against a predicted 6.5, lands at 0.70 — in the
sustainable band. The split-EO populations added on 2026-10-21
(P34, P13, P16) lose 21 – 23 alleles each and still respond in three
different ways. The decomposition runs on all six sites. **Without the diversity trigger from § C.0.a we
would not have known there was a mechanism question to ask.**

**Three competing hypotheses can produce different pollen-
compatibility responses to the same diversity loss.**

1. **Within-class allele spread.** Class I and Class II total
   masses can be identical, but all of Class I's mass may be
   locked in a single allele (drift monoculture), killing
   within-class compatibility (a Class I × Class I cross is
   rejected whenever parents share an allele).
2. **Class I / Class II mass balance.** Raising / lowering the
   Class I share changes how much "between-class rescue" is
   available to Class II mothers (between-class crosses are
   always compatible under the sporophytic model).
3. **Zygosity composition.** The fraction of mothers that carry
   1 / 2 / 3 distinct SRK identities. Multi-identity mothers
   express a bigger set and face more compatible fathers for
   Class II; the effect is more complex for Class I.

#### How the decomposition is built

The same simulation that produced Figure 11 (pollen compatibility
per location) is run **five
times per location**, keeping the observed mother genotypes fixed
and swapping ONE part of the father-drawing distribution at a time
against the species-wide reference. For each scenario, 800
candidate fathers are drawn, sporophytic Class I / II compatibility
is evaluated per mother, averaged over her trials, and then over
all her location's mothers. Bootstrap over 400 resamples of
mothers gives a 95 % CI.

The **only** thing that differs between the five bars for a given
location is the father-drawing distribution:

| Scenario (bar colour) | Father allele frequencies | Father zygosity composition |
|--- |--- |--- |
| **Observed** (black) | Observed local | Observed local |
| **Swap within-class spread** (blue) | **P1 shape within each class, scaled to the observed Class I and Class II totals** | Observed local |
| **Swap Class I / II balance** (red) | **Within-class shape kept observed, class totals rescaled to the species-wide values** | Observed local |
| **Swap zygosity** (yellow) | Observed local | **Species-wide 66 / 32 / 2 % of 1-/2-/3-distinct** |
| **prediction** (green) | Species-wide P1 | Species-wide |

For each single-swap bar, the **driver fraction** = (counterfactual
mean − observed mean) / (species-wide prediction mean − observed
mean). A value near **+100 %** says that one factor alone closes
the entire observed → species-wide gap. A value near **0 %** says
the factor was not involved. A **negative value** means swapping
away from the observed state moves pollen compatibility *away*
from the species-wide prediction — the observed state on that axis
is **buffering** the location.

#### How to read the figure

- **Observed black bar below green**: the location has lost
  pollen compatibility relative to the species-wide prediction.
- **A single-swap bar jumps up toward green**: that factor
  explains most of the gap.
- **A single-swap bar stays next to black**: that factor is not
  involved.
- **A single-swap bar drops below black**: that factor is
  buffering the location — observed state on that axis is more
  favourable than species-wide.

#### Walk-through — the six genotyped sites

| Location | Alleles obs / pred | Pollen compatibility obs → species-wide pred | Within-class spread | Class I / II balance | Zygosity composition |
|---|---|---|---:|---:|---:|
| **EO67 (P6)** | 7 / 6.5 (+0.5) | 0.70 → 0.76 | **+100 %** | −7 % | −13 % |
| **EO70 (P42)** | 6 / 22.5 (−16.5) | 0.54 → 0.70 | **+88 %** | −3 % | −6 % |
| **EO76 (P40 + P41)** | 9 / 29.4 (−20.4) | 0.75 → 0.77 | **+206 %** | +1 % | **−115 %** |
| **P34 (EO18)** | 7 / 29.2 (−22.2) | 0.73 → 0.74 | +172 % | +228 % | **−472 %** |
| **P13 (EO27)** | 5 / 27.8 (−22.8) | 0.36 → 0.69 | +4 % | **+118 %** | +21 % |
| **P16 (EO27)** | 8 / 29.4 (−21.4) | 0.77 → 0.86 | +46 % | +4 % | **+80 %** |

Percentages are shares of each site's gap (species-wide − observed)
closed by that single swap; for P34 the gap is only 0.01, so its
percentages are large but the effects are small. The "observed"
values here use each site's own father zygosity; Figure 13 uses the
species-wide zygosity for fathers, which is why P16 sits above its
prediction there (0.84) and below it here (0.77).

- **P42 — classical drift collapse.** FG024 absorbs 35 % of the
  pool, FG001 another 62 %; four other alleles sit under 1 %.
  A FG024-homozygous Class I mother sees `(1 − 0.35)⁴ ≈ 0.17`
  compatibility, vs ≈ 0.43 at P6 where Class I is split four
  ways. The blue bar jumps up to 0.68 (closing most of the gap);
  red and yellow stay next to black. **Within-class concentration
  explains the deficit; class balance and zygosity are not
  involved.**
- **P6 — looks bad on paper, pollen compatibility survives.**
  Only 7 of 32 alleles present, but those 7 split as **4 Class I
  alleles (FG024 / FG018 / FG023 / FG016, well spread) and 3
  Class II alleles**. Each Class I mother expresses a smaller
  share of Class I mass than at P42, so Class I × Class I
  compatibility stays reasonable. Observed zygosity (32 %
  multi-identity) is slightly higher than species-wide — the
  yellow bar sits just below black, meaning observed zygosity is
  mildly buffering. **Mate limitation at P6 is not as bad as
  the headline 7/32 count suggests.**
- **P40 — predicted near-perfect, observed severely collapsed,
  pollen compatibility still fine.** 9/32 alleles — a collapse
  comparable to P42 — but pollen compatibility drops only 0.03.
  Decomposition exposes the mechanism: **zygosity is actively
  buffering P40** (yellow sits at 0.72, below the black 0.75,
  a −115 % move relative to the gap). The location is 76 %
  homozygous (vs species-wide 66 %), and because its dominant
  Class II allele (FG001 at 46 %) is less dominant than P42's
  FG001 (62 %), those homozygous Class II mothers express a
  small, non-dominant set and face many compatible fathers. Blue
  overshoots the species-wide prediction (within-class diversity
  at P40 is actually better-spread than drift-only P1 would
  deliver).
- **P34 (EO18) — almost no gap here.** With its own zygosity P34's
  compatibility (0.73) matches the species-wide value (0.74):
  within-class spread and class balance would each raise it a
  little, and its own zygosity buffers it (swapping to the
  species-wide zygosity would lower it to 0.68, yellow bar). Its
  shortfall in Figure 13 comes from the comparison with the
  population-level prediction (0.78).
- **P13 (EO27) — a class-balance collapse.** Only 4 % of P13's SRK
  copies are Class I (35 % species-wide): FG001 and FG002 dominate
  and almost every plant expresses Class II, so most crosses are
  between plants that share a dominant Class II allele. Restoring
  the species-wide class balance lifts compatibility from 0.36 to
  0.75 — the whole gap (red bar). This is a different failure from
  P42's within-class monoculture. 12 plants from 7 mothers:
  indicative.
- **P16 (EO27) — zygosity, not allele loss.** 53 % of P16's plants
  express three or more distinct alleles (about 2 % species-wide).
  Fathers that express more identities match more mothers, so the
  site's own zygosity pulls compatibility down to 0.77; with
  species-wide zygosity it reaches 0.84 (yellow bar). Many of these
  plants are seedlings, whose genotypes combine two parents. 15
  plants from 7 mothers: indicative.

#### Methodological take-home

A location's SRK diversity gap tells us drift has happened; the
pollen-compatibility response is **not** a monotonic function of
how many alleles were lost. It depends on **which** alleles were
lost, **how dominant** the remaining alleles are, and **how multi-
identity** the mothers are. Capturing this requires an accurate
sporophytic Class I / II model plus empirical zygosity (§ A.7–A.8)
— a simpler diploid gametophytic approximation would collapse all
three channels into a single "effective diversity" number and
mis-predict P40 outright. The hypothesis decomposition therefore
doubles as validation of the model's mechanistic structure: it
successfully resolves cases where the three channels pull in
different directions.

#### Biological take-home — the deme can buffer intense drift

The three single-swap counterfactuals (within-class spread, Class I /
Class II balance, tetraploid zygosity) are not independent model
knobs — they are three **buffering mechanisms** stacked on top of
each other that let a deme keep breeding even after severe allele
loss:

1. **Class I dominance (26 of 32 Fgs).** Every between-class cross
   is compatible by construction, and Class I dominates Class II
   within a plant. A mother carrying any Class I allele has her
   compatibility set primarily by that allele, so losing *rare*
   Class I alleles barely moves the mean (the common Class I
   alleles still carry the pool).
2. **Frequency-dependent selection is forgiving when drift
   collapses onto *common* alleles.** The classical SI catastrophe
   (Lawrence 2000; Castric & Vekemans 2004) is a deme that collapses
   onto so few alleles that most mothers share them all → crash.
   But at P42 the 6 surviving Fgs are the species's *common* ones
   (FG001 at 41 %, FG024 at 35 % locally) — the opposite of the
   worst case. Different mothers still carry different combinations,
   so the pollen pool still finds compatible targets.
3. **Tetraploid zygosity actively buffers at homozygous-rich sites.**
   At P40, 76 % of mothers are homozygous: a homozygous mother
   expresses a single Fg at her stigma, which pollen fathers can
   more easily avoid matching than a heterozygous mother's larger
   expressed set. The yellow bar at P40 in Figure 14 **drops below
   the observed black bar** — the quantitative signature of this
   mechanism.

**What this reverses in the standard conservation-genetics
intuition.** The usual story — *lose SRK alleles → mating failure*
— treats allele count as the bottom line. In a **tetraploid +
sporophytic Class I / II system with the specific class imbalance
LEPA has**, the mating-level consequence of allele loss is
**actively buffered** until the pool either (i) collapses onto a
single class or (ii) homozygosity becomes so extreme that the
expressed-set channel constrains compatibility. The 32 → 6 drift
collapse at P42 is severe by any count, and yet the pollen pool
still works at 70 % compatibility (close to the species mean 0.78).
**The system has structural redundancy that diversity-counting
alone cannot see.** The next-sampling-campaign data at the B1 pilot
pair will tell us whether any location is near either tipping point.

**Workflow placement.** This analysis is a follow-through on § C.0.a's
diversity trigger, not a stand-alone pollen-compatibility test.
When a location's observed allele diversity sits noticeably below
the prediction, run `step26_partC_hypothesis_decomposition.py`
to see which of the three channels is driving the downstream
pollen-compatibility response.

**Outputs.**

- [`step30_B_partC_hypothesis_decomposition.tsv`](tables/step30_B_partC_hypothesis_decomposition.tsv)
  — one row per clean-overlap location with the diversity trigger
  (observed vs predicted distinct allele counts), observed +
  three single-swap counterfactual pollen compatibility estimates
  with bootstrap 95 % CIs, observed Class I mass, observed zygosity
  distribution, and the three "driver fraction" columns.
- [`step30_B_partC_hypothesis_decomposition.png`](figures/step30_B_partC_hypothesis_decomposition.png)
  — grouped bar chart per location (see caption below).

<a id="fig-14"></a>
![Figure 14](figures/step30_B_partC_hypothesis_decomposition.png)

**Figure 14.** Competing-hypothesis decomposition of pollen compatibility at the six genotyped sites — P6 (= EO67), P42 (= EO70), P40 + P41 (= EO76), P34 (EO18), P13 and P16 (EO27). For each site the figure shows five pollen-compatibility values from the same simulation (sporophytic Class I / II + empirical zygosity, 800 candidate fathers per mother, observed mother genotypes held fixed), each differing in which part of the father-drawing distribution is swapped to the species-wide reference. **Green — predicted:** everything species-wide. **Grey — observed:** fathers drawn from the site's allele frequencies and zygosity. **Blue — swap within-class spread:** keep the observed Class I and Class II totals, reshape the spread within each class to the species-wide pattern (the "drift monoculture" channel). **Red — swap Class I / II balance:** keep the within-class shape, rescale the two class totals to species-wide (between-class rescue). **Yellow — swap zygosity:** keep the allele frequencies, swap father zygosity to species-wide 66 / 32 / 2 % (the expressed-set channel). Star = the swap that moves compatibility furthest, the dominant explanation of the gap. Error bars = bootstrap 95 % interval (400 mother resamples). Traffic-light bands in the background; dashed grey line = species-wide pollen compatibility 0.78. **Reading rules:** a coloured bar rising toward green means that factor causes the gap; staying next to grey means it is not involved; dropping below grey means it is **buffering** the site. **P42** — within-class monoculture (FG024 at 35 % of the pool). **P6** — small deficit, within-class spread explains it. **EO76** — small deficit despite 9 of 32 alleles; its zygosity buffers it. **P34** — almost no gap; zygosity buffering. **P13** — class balance: only 4 % of its SRK copies are Class I, and restoring the species-wide balance closes the whole gap. **P16** — zygosity: 53 % of its plants express ≥ 3 identities, which lowers compatibility. P13 and P16 rest on 12 and 15 plants. Source: `step26_partC_hypothesis_decomposition.py`.

**Note — why P1 is empirical, not uniform.** A uniform-frequency P1
(Dirichlet(α = 1) over 32 Fgs, equivalent to the P0 "uninformative"
prior kept as a reference baseline) would represent a **neutral null
with no drift history**. That is not LEPA: decades of habitat loss
have already pushed the species through drift, so empirical P1 — the
observed Fg frequencies in the Canu-amplicon 263-individual inventory
— is the realistic *starting point* against which per-location
further drift is measured. Panel C makes this concrete: P42's local
pool has moved *past* P1 toward FG001 dominance, so even the already-
skewed empirical prior underestimates how collapsed the local pool
is. A uniform P1 would start from a flat distribution and declare
every observed location "drifted", which is both less informative
(the species-wide signal is real) and less actionable (the baseline
would not reflect how LEPA actually enters the modelling frame).

### C.0.c Hypothesis test — can the diversity gap be closed by tightening the deme radius?

**Why this test exists.** Figure 13 shows the diversity prediction
over-estimates the observed SRK allele count at every genotyped site
(predicted 22 – 29 vs observed 5 – 9 at five of the six sites). § C.0.b above decomposes *why
pollen compatibility nevertheless tracks observation* at some sites
and not others; here we ask whether the diversity gap at EO70 (P42)
and EO76 (P40 + P41) could be an artefact of the operational deme
being too wide at 50 m. The sweep re-partitions each location's
slick spots, so it runs on these two EO sites (P6 is invariant
across radii); the split-EO populations P34, P13 and P16 are not
swept. Two competing scenarios:

- **H1 — the 50 m operational deme is too wide.** Realised gene
  flow is tighter; a smaller radius should partition each location
  into more small demes with (within the current model) independent
  drift histories drawn from the species-wide prior P1. If H1 is
  right, the gap should shrink as radius tightens.
- **H2 — per-deme drift history exceeds the species-wide prior
  P1.** Each deme has drifted *past* what P1 captures, so no
  spatial repartition at any radius closes the gap under the
  current model. The gap should stay flat across the sweep.

**Scope — what the test can and cannot do.** The sweep rebuilds
the deme partition at radii 10, 25, 50, 75, 100, 150 m and re-runs
the per-deme diversity simulation (2000 replicates per combination).
Each sub-deme draws alleles from P1 — i.e. the test varies the
*spatial partition* while holding the *drift model* fixed. **It
can only reject H1 within the P1 drift assumption; it does not
simulate H2 directly** (that would require per-deme empirical
frequency vectors we do not have yet — the next seed-genotyping
campaign is designed to measure them).

**P6 is omitted from the figure.** It has only 2 events in the
2025 record, placed > 150 m apart, so the deme partition is
invariant across every tested radius (2 demes everywhere). The
sweep is non-informative there, and P6's observed count is
inside the 95 % CI of the prediction already — no gap to
explain. The figure restricts to **EO70 (P42)** and **EO76 (P40 +
P41)**, which do show meaningful fragmentation across the sweep (EO70:
1 → 3 demes; EO76: 4 → 16 demes as radius tightens).

**Result — H1 is rejected within the P1 drift assumption at EO70
and EO76** ([Figure 15](#fig-15)).

| Location | Observed | Predicted at r = 10 m | at 50 m | at 150 m | Deme count 10 m / 50 m / 150 m |
|---|:---:|:---:|:---:|:---:|:---:|
| **EO70 (P42)** | 6 | 22.5 [19, 27] | 22.6 [19, 26] | 22.5 [19, 26] | 3 / 1 / 1 |
| **EO76 (P40 + P41)** | 9 | 29.4 [26, 32] | 29.4 [27, 32] | 29.4 [26, 32] | 16 / 4 / 4 |

The prediction **does not move** across the sweep: even fragmenting
EO76 into 16 small demes at 10 m, or P42 into 3 demes at 10 m,
leaves the union-of-per-deme pools at essentially the same 22 – 29
Fgs. The ~17-allele gap at P42 and ~20-allele gap at P40 do not
close. Full table at
[`step30f_srk_diversity_radius_sweep.tsv`](tables/step30f_srk_diversity_radius_sweep.tsv).

**Why the sweep is flat — the mechanics.** The location-level
prediction is the **union** of per-deme Fg sets, not a sum. When
a radius change re-partitions a location into more or fewer
sub-demes, the **total number of allele copies sampled stays
fixed at 4 × N_fertile_total** — the sweep only redistributes
those draws across sub-demes. Each sub-deme on its own may miss
rare Fgs under P1 (where FG001 ≈ 41 %, five more Fgs ≈ 25 %
combined, 20+ rare Fgs each < 2 %), but a rare Fg at frequency
*p* now has multiple sub-demes to appear in, and the union captures
it with probability `1 − (1 − p)^(4 × N_total)` — governed by the
location **total**, not the partition. Concretely at P42's 10 m
extreme the 161 plants split into sub-demes of 78 / 65 / 18 that
each recover ~26 / 25 / 19 Fgs independently, but **their union
still delivers 28**, same as the single 161-plant deme at 50 m.
Spatial sub-partitioning washes out at the location level once the
total pool saturates P1 (above ~50–100 adults in LEPA). The
partition would only matter at a location whose total `N_fertile`
sits below that saturation ceiling and whose sub-deme sizes drop
into the low single digits — neither applies at P42 or P40.
**Under P1 drift, varying the deme partition therefore cannot
close the diversity gap at any radius** — only a drift model that
differs per deme can (which is what the seed-genotyping campaign
will measure).

**Interpretation — what we can and cannot say.**

1. **The 50 m operational deme is defensible as a first-pass
   partition.** No tighter radius would reconcile prediction with
   observation under the current model.
2. **The gap cannot be a pure spatial-partitioning artefact.** The
   only way it closes is if drift histories *inside* each deme
   have diverged from P1, i.e. per-deme empirical priors — that
   is the H2 scenario, which requires seed-genotyping to test.
   **Biologically: the sweep rules out a geometric explanation
   and leaves the within-deme biological one — prolonged drift,
   founder effects, or local bottlenecks that have stripped
   alleles beyond what the species-wide prior encodes. The deme
   is still the right unit of inference; what the data say is
   that the drift history inside each deme is more severe than a
   single species-wide frequency vector captures, which is the
   central prediction of the fragmentation → drift → mate-
   limitation causal chain this framework was built to test.**
3. **The direction of error is favourable.** Over-predicting at a
   radius that is already a geographic upper bound says the real
   demes are **at most** 50 m and could be tighter. A future
   refinement can only add demes (and therefore mothers per
   location), giving us finer-scale information for free. The
   opposite — under-predicting ⇒ deme too narrow ⇒ false-positive
   fragmentation — would waste sampling on sub-demes that do not
   exist.

<a id="fig-15"></a>
![Figure 15](figures/step30f_srk_diversity_radius_sweep.png)

**Figure 15.** SRK diversity gap (y = predicted − observed distinct SRK alleles) vs pollinator radius at **EO70 (P42)** and **EO76 (P40 + P41)**. x = radius used to rebuild the deme partition (10 → 150 m, log scale). Lines + 95 % CI bands = per-location per-radius simulation (2000 replicates per combination). Right-side labels show the gap at the largest radius; parenthetical observed counts come from the Phase 4 adult SRK genotypes (a single measurement, not a sweep output — shown as labels, not horizontal lines). Zero line = perfect match. Vertical dotted line at 50 m marks the current operational deme. H1 (deme too wide) would predict the gap to shrink toward zero as radius tightens; it stays flat at +22 and +23 across all radii, including the 10 m extreme that fragments P42 into 3 demes and EO76 into 16 — rejecting H1 within the P1 drift assumption. P6 (= EO67) is omitted (deme partition invariant across the sweep, test non-informative). Source: `step27_radius_sweep_H1.py`.


## Phase B — Closing the causal chain with seed data

**Test 1 — Mate-limitation regression (§ C.1 of the long doc).**
Once seed genotypes are back, per-location observed P_compat is
computed from real mother-father pairings. The regression
`per-mother seed set ~ β₁ · P_compat + β₂ · mating_neighbourhood + …`
decomposes the fragmentation effect into two pathways:

- **β₁ (P_compat effect)** measures the strength of the **complete
  drift → mate-limitation branch** — the reproductive consequence of
  the fragmentation-driven Fg loss and skew that Phase A predicted.
- **β₂ (mating-neighbourhood effect)** measures whether fragmentation
  has **direct effects on seed set that are NOT mediated through Fg
  composition** (e.g. fewer neighbouring adults means fewer pollinator
  visits, independent of which alleles they carry).

Together the two coefficients tell us how much of fragmentation's
reproductive cost flows through drift, and how much bypasses it. A
large β₁ + small β₂ = the drift chain fully explains the fragmentation
effect; a large β₂ + small β₁ = fragmentation reduces seed set through
pollinator-behaviour pathways we haven't modelled; both large = both
channels contribute. All three interpretations can be pre-registered
before seed data arrive.

---

## Recommended field pilot — BL3 within-BL stable pair (P36 + P34)

Before running Phase B across all 44 populations, we recommend a
**two-population pilot within Bottleneck Lineage 3 (BL3)** — one
**SMALL stable** population (**P36 = EO25-B_21**, 26 adults
across both years, 7 + 8 = **15 mothers in DB**) paired with one
**BIG stable** population (**P34 = EO18-7 + EO18-8**, 1 673 adults
across five locationIDs, 115 + 221 = **336 mothers in DB**).
Both pass the stability test (both-year occupancy, overlapping
95 % 2025 vs 2026 CIs on diversity and pollen compatibility,
`N_fert_total ≥ 20`). Updated 2026-10-07 from the earlier
BL5-fragmentation-contrast design; see § Across-year population
dynamic for the stable-population inventory and the Executive
summary **Next test** block for the full D1/D2/D3 option table.

**Why BL3.** Across the three BLs that hold ≥ 2 stable populations
(BL1, BL2, BL3), **BL3 holds the single largest well-stocked
population** (P34, 1 673 adults with 336 mothers with cleaned seed
in the LEPA DB) **and** a stable SMALL population at the other end
of the size spectrum (P36, 26 adults). The pair-level mother budget
of **351 mothers across 2025 + 2026** is the largest available
within-BL stable contrast in the dataset and roughly **1.8 × the
next-best option** (D3 at 190 mothers, BL1 P6 + P16). Holding BL identity constant means observed
deviations are a within-BL small-vs-big signal, not a between-BL
drift-history confound.

**What this specific pair tests.**

- **P34 EO18-7 + EO18-8 — large stable regime** (1 673 adults
  spread across 5 locationIDs and several demes). The model
  predicts a sustainable population at the species mean; this
  side is the cleanest test of the model's "ample pool →
  compatible mating" prediction. **Caveat:** "stable" here means
  overlapping predicted diversity and compatibility between years;
  the census is not constant (523 → 1 150 adults, 16 → 21 demes),
  and seed yield was above expectation in 2025 (fold 1.28; at
  expectation once the 2026 year effect is allowed for, so the call
  is year-sensitive). No large population in the dataset has both a
  constant census and seed yield at expectation in both years; P34
  is kept for the within-BL design and its 336 mothers.
- **P36 EO25-B_21 — small stable regime** (26 adults). Small
  enough that drift should erode the pool measurably per deme,
  yet stable enough to carry a mother budget. If the
  per-deme simulation predicts a materially lower pollen
  compatibility at P36 than P34, and the observed seed-set data
  track that gap, the within-BL small-vs-big signal is confirmed.

**Sampling the pair.** P36 holds 15 mothers with seed across the two
years and P34 336. The two-step design samples 7 of P36's and 84 of
P34's, in the demes chosen to recover 90 % of each population's SRK
alleles — 91 mothers, ≈ 1 360 seedlings; see [`Sampling_design.md`](Sampling_design.md) § 4.

**Pair evaluations and stable-population inventory**:
[`step30g_populations_classified.tsv`](tables/step30g_populations_classified.tsv)
is the authoritative source of the trend classification plus the
`populationID → locationCode` map; the three within-BL stable
pair options (D1 / D2 / D3) are summarised in the Executive
summary.

---

## Deliverables and where to look

| What | File | Notes |
|---|---|---|
| P1 prior over 32 Fgs | [`tables/step30_A_prediction_prior_frequencies.tsv`](tables/step30_A_prediction_prior_frequencies.tsv) | Species-wide Fg frequencies + Dirichlet α |
| Predicted per-location diversity | [`tables/step30_A_prediction_location_diversity.tsv`](tables/step30_A_prediction_location_diversity.tsv) | Species-wide + local coverage, 95 % CI |
| Predicted per-50 m-component diversity | [`tables/step30_A_prediction_component_diversity.tsv`](tables/step30_A_prediction_component_diversity.tsv) | Pool size, allocated mothers/seeds, coverage per component |
| Predicted per-location P_compat | [`tables/step30_A_prediction_location_pcompat.tsv`](tables/step30_A_prediction_location_pcompat.tsv) | Sporophytic + empirical zygosity |
| Predicted per-50 m-component P_compat | [`tables/step30_A_prediction_component_pcompat.tsv`](tables/step30_A_prediction_component_pcompat.tsv) | Exposes struggling sub-components inside a sustainable-mean location |
| Fragmentation indices | [`tables/step30_A_fragmentation_per_event.tsv`](tables/step30_A_fragmentation_per_event.tsv), [`_per_population.tsv`](tables/step30_A_fragmentation_per_population.tsv) | Pure spatial, no allele frequencies |
| Event → population → deme | [`tables/step29a_event_to_deme.tsv`](tables/step29a_event_to_deme.tsv) | Canonical "which events share a pollen pool?" |
| **Sampling design, sampling check and lab lists** | [`Sampling_design.md`](Sampling_design.md) | Two-step plan; germplasmIDs per deme for the model populations |
| Seed yield vs predictions | [`tables/step30i_seed_yield_vs_predictions.tsv`](tables/step30i_seed_yield_vs_predictions.tsv) + [`_tests.tsv`](tables/step30i_seed_yield_vs_predictions_tests.tsv) | Fold of expectation next to predicted compatibility and diversity; Spearman tests |
| Background EO-scale diagnostic for split EOs (§ C.0) | [`tables/step30_C_pcompat_validation_at_eo.tsv`](tables/step30_C_pcompat_validation_at_eo.tsv) | EO18, EO25, EO27 (the three EOs whose location-scale remap is still pending); no figure |
| **Part C anchor at location scale (§ C.0.a)** | [`tables/step30_B_partC_clean_overlap_per_location.tsv`](tables/step30_B_partC_clean_overlap_per_location.tsv) + [`_fg_frequencies.tsv`](tables/step30_B_partC_clean_overlap_fg_frequencies.tsv) | **Six genotyped sites — EO67 (= P6), EO70 (= P42), EO76 (= P40 + P41) and P34 (EO18), P13, P16 (EO27); 246 plants; observed vs predicted pollen compatibility + SRK diversity + no-drift upper bound** |
| **Hypothesis decomposition (§ C.0.b)** | [`tables/step30_B_partC_hypothesis_decomposition.tsv`](tables/step30_B_partC_hypothesis_decomposition.tsv) | **Decomposes the pollen-compatibility deviation at each clean-overlap location into within-class spread / Class I-II balance / zygosity contributions; triggered by the SRK diversity gap** |

**Key figures** (ordered by appearance in this doc).

- **Figure 1** — [`step30_A_radius_sensitivity.png`](figures/step30_A_radius_sensitivity.png) — why 50 m is the right primary pollinator radius.
- **Figure 2** — [`step29a_populations_overview.png`](figures/step29a_populations_overview.png) — the 44 populations across the Snake River Plain, colour-coded by 2025 + 2026 occupancy.
- **Figure 3** — [`step30h_dendrogram.png`](figures/step30h_dendrogram.png) — Ward's D2 dendrogram at the population level; the *k* = 5 cut defines the five Bottleneck Lineages.
- **Figure 4** — [`step30h_silhouette_curve.png`](figures/step30h_silhouette_curve.png) — silhouette curve; optimum *k* = 8, *k* = 5 kept for continuity with the EO-level lineages (0.03 lower).
- **Figure 5** — [`step30h_overview_by_BL.png`](figures/step30h_overview_by_BL.png) — Snake River Plain overview with populations and convex hulls coloured by BL1 … BL5.
- **Figure 6** — [`step30g_populations_classified.png`](figures/step30g_populations_classified.png) — the 44 populations grouped by across-year trend class (crash / stable / growth / ambiguous / 2025_only / 2026_only), sorted by BL within each class.
- **Figure 7** — [`step30g_across_year_scatter.png`](figures/step30g_across_year_scatter.png) — predicted SRK diversity and pollen compatibility, 2025 vs 2026, per both-year population.
- **Figure 8** — [`step29a_demes_per_population_year.png`](figures/step29a_demes_per_population_year.png) — per-population per-year deme-size distribution (log-scale deme size with N = 1 / 8 / 32 reference lines, 2025 open vs 2026 filled).
- **Figure 9** — [`step30h_pred_diversity.png`](figures/step30h_pred_diversity.png) — predicted SRK allele diversity per population, 2025 (open) vs 2026 (filled) overlaid, grouped by new BL.
- **Figure 10** — [`step30_A_si_model_schematic.png`](figures/step30_A_si_model_schematic.png) — sporophytic Class I / Class II dominance in tetraploid LEPA: dominance within one plant, between-plant recognition, and the compatibility rule by cross type.
- **Figure 11** — [`step30h_pred_pcompat.png`](figures/step30h_pred_pcompat.png) — predicted pollen compatibility per population, 2025 (open) vs 2026 (filled) overlaid, traffic-light bands, grouped by new BL.
- **Figure 12a** — [`step30i_species_calibration.png`](figures/step30i_species_calibration.png) — § 30.3 across-species build: calibration scatter (observed vs expected per plant, coloured by BL) + predictor-selection CV table + expectation-model summary.
- **Figure 12b** — [`step30i_size_seed_population.png`](figures/step30i_size_seed_population.png) — § 30.3 within-population test: observed seed yield per population vs size expectation, 2025 (open) vs 2026 (filled), all 44 populations with ≥ 2 plants; **red halo** = FDR-sig below (mate-limited), **green halo** = FDR-sig above (buffered surplus) for that (population, year). Priority list in [`step30i_conservation_priority.tsv`](tables/step30i_conservation_priority.tsv) (one row per populationID; one restoration status built from both years, ranked by `priority_rank`).
- **Figure 12c** — [`step30i_seed_yield_vs_predictions.png`](figures/step30i_seed_yield_vs_predictions.png) — seed yield against predicted pollen compatibility and SRK diversity.
- **Figure 13** — [`step30_B_partC_clean_overlap.png`](figures/step30_B_partC_clean_overlap.png) — § C.0.a Part C anchor at P6 (EO67) / P42 (EO70) / P40 (EO76) (SRK diversity + pollen compatibility + per-allele drift residual).
- **Figure 14** — [`step30_B_partC_hypothesis_decomposition.png`](figures/step30_B_partC_hypothesis_decomposition.png) — § C.0.b competing-hypothesis decomposition.
- **Figure 15** — [`step30f_srk_diversity_radius_sweep.png`](figures/step30f_srk_diversity_radius_sweep.png) — § C.0.c SRK-diversity radius sweep (deme-size vs residual-drift test).
- Candidate-population zoom panels — [`step29a_candidate_populations_zoom.png`](figures/step29a_candidate_populations_zoom.png) + [`step29a_crash_populations_zoom.png`](figures/step29a_crash_populations_zoom.png) (supporting plots for the LARGE + SMALL + crash candidates).

---

## For technical detail

Formulas, code, output schemas, and design decisions are documented
in the archived long doc [`archive/Sampling_and_prediction.md`](archive/Sampling_and_prediction.md) (no longer maintained; this outline takes precedence where they differ):

- § A.3 — species-wide prior + two-generation model.
- § A.5 — fragmentation of pollen flow, event vs location scales.
- § A.7 — sporophytic Class I / Class II SI biology.
- § A.8 — finite-population P_compat formulas (Case A / Case B).
- § B.4.2 — fragmentation-aware sampling derivation.
- § C.0 — EO-level empirical validation (the § C.0 above is the
  compact version of this section).
- § C.1 — the mate-limitation regression (Phase B).

---

## References

The population-genetics scaffolding for the operational-deme
framework and the fragmentation → drift → mate-limitation causal
chain. BibTeX entries live in [`SRK_mate_limitation_references.bib`](SRK_mate_limitation_references.bib).

- **Aguilar, R., Quesada, M., Ashworth, L., Herrerías-Diego, Y. & Lobo, J. (2008).** Genetic consequences of habitat fragmentation in plant populations: susceptible signals in plant traits and methodological approaches. *Molecular Ecology* 17, 5177–5188.
- **Castric, V. & Vekemans, X. (2004).** Plant self-incompatibility in natural populations: a critical assessment of recent theoretical and empirical advances. *Molecular Ecology* 13, 2873–2889.
- **Freckleton, R.P. & Watkinson, A.R. (2002).** Large-scale spatial dynamics of plants: metapopulations, regional ensembles and patchy populations. *Journal of Ecology* 90, 419–434.
- **Hanski, I. (1998).** Metapopulation dynamics. *Nature* 396, 41–49.
- **Hardy, O.J. & Vekemans, X. (1999).** Isolation by distance in a continuous population: reconciliation between spatial autocorrelation analysis and population genetics models. *Heredity* 83, 145–154.
- **Levins, R. (1969).** Some demographic and genetic consequences of environmental heterogeneity for biological control. *Bulletin of the Entomological Society of America* 15, 237–240.
- **Honnay, O. & Jacquemyn, H. (2007).** Susceptibility of common and rare plant species to the genetic consequences of habitat fragmentation. *Conservation Biology* 21, 823–831.
- **Lawrence, M.J. (2000).** Population genetics of the homomorphic self-incompatibility polymorphisms in flowering plants. *Annals of Botany* 85 (Suppl. A), 221–226.
- **Levin, D.A. & Kerster, H.W. (1974).** Gene flow in seed plants. *Evolutionary Biology* 7, 139–220.
- **Schierup, M.H., Vekemans, X. & Christiansen, F.B. (1998).** Allelic genealogies in sporophytic self-incompatibility systems in plants. *Genetics* 150, 1187–1198.
- **Vekemans, X. & Hardy, O.J. (2004).** New insights from fine-scale spatial genetic structure analyses in plant populations. *Molecular Ecology* 13, 921–935.
- **Wright, S. (1943).** Isolation by distance. *Genetics* 28, 114–138.
- **Wright, S. (1946).** Isolation by distance under diverse systems of mating. *Genetics* 31, 39–59.
