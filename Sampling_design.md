# Sampling design — SRK diversity of the D1 model populations

**What this document is.** The agreed plan for measuring **SRK allele
diversity** in the two D1 model populations by genotyping seedlings and
their mothers: the goal (§ 1), how the populations were chosen (§ 2),
the sampling strategy at a glance (§ 3), the research behind the design
(§ 4), and the sampling itself — the D1 plan, deme maps and lab sheet
(§ 5). Sampling to measure pollen compatibility (mate limitation) is a
separate step, designed next (§ 6). The wider project — populations, Bottleneck Lineages (BLs),
predictions and the seed-yield test — is described in
[`Project_outline.md`](Project_outline.md). An earlier, more detailed
version of this document is kept, no longer maintained, in
[`archive/Sampling_design_detailed.md`](archive/Sampling_design_detailed.md).

**For the lab team:** the mothers to germinate are listed in
**[`tables/step30_D_model_design_lab_sheet.tsv`](tables/step30_D_model_design_lab_sheet.tsv)**
(§ 5.2).

---

## 1. Goal

**Estimate SRK allele diversity at population level, and test the
hypothesis that smaller populations hold less SRK diversity.**

LEPA is self-incompatible: a plant can only be fertilised by a mate
carrying different SRK alleles. The project follows one causal chain —
**fragmentation → genetic drift → mate limitation or SI escape**
([`Project_outline.md`](Project_outline.md)). Its second link is the
one tested here: genetic drift removes SRK alleles faster in small
populations, so **a small population should hold fewer SRK alleles
than a big one with the same history.**

**Why it needs testing.** The genotypes we already hold show far less
SRK diversity than the project's model expects: at five of the six
genotyped sites we find 5 – 9 alleles where the model, which treats each
population as a random sample of the species, predicts 22 – 29. Fitted to
those six sites, the model suggests that populations behave as if
descended from about a dozen founder plants *whatever their size today*
— which would mean current size predicts SRK diversity poorly. The six
sites differ in lineage and history, though, so they cannot separate
size from history. A small and a big population **of the same BL** can.

**What the test looks for.**
- **P36 (small) holds clearly fewer SRK alleles than P34 (big)** →
  supports the hypothesis: current population size limits SRK
  diversity.
- **P36 holds about as many as P34** → SRK diversity was set by past
  bottlenecks shared by the lineage, not by current size.

**What the sampling must deliver.** For each population and each season
(2025 and 2026 seed), seedlings that recover **90 % of the population's
SRK alleles**, so that a difference between P36 and P34 reflects the
populations, not unequal sampling effort — with as few mothers as
possible and spread across each population.

---

## 2. How the model populations were selected

**Approach — a small and a big population in the same BL.** A BL is a
group of populations with a shared fragmentation history. Comparing a
small and a big population within one BL holds that history constant,
so a difference in SRK diversity between them is a size effect, not a
difference in lineage.

**Criteria.** Both populations must:
1. belong to the same BL;
2. be present in both seasons (2025 and 2026), with at least 20 fertile
   plants over the two years;
3. be **stable** between seasons — the 95 % intervals of their
   predicted SRK diversity and of their predicted pollen compatibility
   overlap between 2025 and 2026 (a population that crashed or boomed
   between seasons would not);
4. contrast strongly in size;
5. hold enough collected mothers to sample every deme they need, in
   both seasons.

**Candidates.** Three BLs hold a qualifying pair; BL4 and BL5 do not.

| Pair | BL | Small population (fertile plants, 2025 + 2026; mothers collected) | Big population (same) | Size ratio | Status |
|---|:-:|---|---|:-:|---|
| **D1** | BL3 | **P36** (EO25-B) · 26 plants · 15 mothers | **P34** (EO18-7 + EO18-8) · 1 673 plants · 354 mothers | **64×** | **Both stable; largest mother pool** |
| D2 | BL2 | P21 (EO26-3 + EO26) · 191 plants · 70 mothers | P19 (EO8) · 876 plants · 215 mothers | 5× | Neither stable (P21 ambiguous, P19 growing); weakest contrast |
| D3 | BL1 | P6 (EO67) · 22 plants · 10 mothers | P16 (EO27-1 + EO27RT) · 1 991 plants · 306 mothers | 90× | Both stable; small side has 10 mothers |

**Why D1 is best for the goal.**
- **A sharp size contrast** (64×) within one BL.
- **Both populations are stable** across the two seasons.
- **The most material:** 369 mothers collected over the two seasons
  (D3: 316, D2: 285), including 15 on the small side. Every deme of both
  populations can be sampled to its target in both seasons.
- **A head start:** 32 plants of P34 are already genotyped (7 SRK
  alleles found), so the 2025 – 2026 seedlings extend an existing
  dataset.
- **D3 is a close second** and the natural follow-up: a slightly larger
  size contrast, both sides stable and both already partly genotyped
  (P6: 37 plants, P16: 15), but fewer mothers on the small side (10)
  for later work on compatibility. D2's contrast is weak and neither
  side is stable.

**Caveat.** P34's census is not constant (523 → 1 150 fertile plants,
16 → 21 demes) — no large population in the dataset is both constant
in census and stable in its predictions. Sampling both seasons
separately keeps that change visible.

---

## 3. The sampling strategy at a glance

**Key criteria**

| Criterion | Value |
|---|---|
| Unit of sampling | the **deme** — slick spots within 50 m of each other, sharing pollinators |
| Goal per population and season | recover **90 %** of the population's SRK alleles |
| Seeds per mother | **25** germinated → about **15 seedlings** genotyped |
| Mothers per deme | scaled to deme size (1 for 2 – 5 plants, 2 for 10 – 20, 5 for 50, 9 for 100, 14 for 200, 23 for 400), **at least one per slick spot** |
| Demes left out | demes with a single fertile plant |
| Spacing between sampled demes | at least **100 m** |
| Choice of mothers in a deme | **at random** among cleaned seed lots of at least 25 seeds |

**Steps**
1. **Size each deme's sample** — mothers needed to recover 90 % of the
   deme's SRK alleles (§ 4.1).
2. **Choose demes across the population** — demes at least 100 m apart,
   the most new alleles per mother first, until 90 % of the
   population's alleles are recovered; in small P36, all demes (§ 4.2).
3. **Draw the mothers at random** within each chosen deme (§ 4.3).
4. **Germinate 25 seeds per mother and genotype about 15 seedlings**,
   recording each seedling against its mother.

**The result:** **91 mothers, 2 260 seeds, about 1 356 seedlings**
over the two seasons, all from cleaned seed lots ready to germinate
(§ 5.1).
- Maps of the sampled demes: **[Figure 2](#fig-2)** (P36, small,
  [`figures/step30_D_model_design_demes_schematic_P36.png`](figures/step30_D_model_design_demes_schematic_P36.png))
  and **[Figure 3](#fig-3)** (P34, big,
  [`figures/step30_D_model_design_demes_schematic_P34.png`](figures/step30_D_model_design_demes_schematic_P34.png)).
- The list of mothers for the lab:
  **[`tables/step30_D_model_design_lab_sheet.tsv`](tables/step30_D_model_design_lab_sheet.tsv)**
  (§ 5.2).

---

## 4. The research behind the design

### 4.1 Mothers per deme

A seedling carries 2 SRK alleles from its mother and 2 from its father,
a plant of the mother's deme, so a mother's seedlings reveal her own
alleles and alleles of the plants around her. Seeds germinate at about
60 %, so 25 seeds give about 15 seedlings; beyond 15 seedlings a mother
adds little.

Large demes hold more alleles, including rare ones, and need more
mothers. The number is set by simulation. Each deme is simulated plant
by plant:
- tetraploid plants (4 SRK copies), most carrying a single allele in all
  four copies, as observed in LEPA;
- alleles drawn from the species-wide frequencies, so small demes lose
  rare alleles by chance;
- fathers that must be compatible with the mother under LEPA's
  self-incompatibility system.

The target is the smallest number of mothers whose seedlings recover,
on average, 90 % of the alleles the deme's plants carry:

| Fertile plants in the deme | 2 – 5 | 10 – 20 | 50 | 100 | 200 | 400 |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| Mothers needed | 1 | 2 | 5 | 9 | 14 | 23 |

Every slick spot of a sampled deme also gets at least one mother, since
an isolated slick spot may carry alleles found nowhere else. A deme with
a single fertile plant is left out: it has no mate in its deme.

### 4.2 Choosing demes across the population

Sampling every deme of P36 and P34 would take 158 mothers. Instead,
demes are added one at a time, the one adding the most new SRK alleles
per mother first, until the selected demes are expected to recover 90 %
of the whole population's alleles. Sampled demes must be **at least
100 m apart**, so each is a separate pollen pool and the sample spreads
across the population. Because adding the cheapest deme first can, through
the spacing rule, block a deme the population needs later, the search
is repeated from every possible starting deme and the plan with the
fewest mothers is kept. In small P36 both demes are taken.

### 4.3 Mothers chosen at random

Within a chosen deme, mothers are drawn at random among those whose
cleaned seed lot holds at least 25 seeds — one per slick spot first,
then the rest of the deme's number. Mothers are deliberately **not**
chosen by how much seed they set: a mother sets much seed when she finds
many compatible mates, so choosing on seed would favour some allele
combinations and misrepresent the population's SRK diversity. The draw
is reproducible (a fixed random seed per population, season and deme).

### 4.4 One adjustment in 2025

P34's largest 2025 deme (deme 5, 158 plants, a third of the population
that season) holds 9 collected mothers where § 4.1 asks for 12. Leaving
it out would leave a third of the population unrepresented and the plan
short of 90 %, so it is sampled with all 9 of its mothers; its own allele
estimate is a little less precise than the other demes'. P34 was not
sampled less intensively overall (22 mothers collected per 100 plants in
2025, 21 in 2026); this one deme was. 2026 needs no adjustment.

### 4.5 Does the plan recover the SRK diversity? — simulations

The plan was tested on simulated populations built from the real deme
sizes of P36 and P34 ([Figure 1](#fig-1)):

- **Each sampled deme recovers most of its alleles** (Figure 1 B): the
  alleles its mothers' seedlings reveal (red) sit close to the alleles
  the deme holds (grey).
- **Each population reaches the 90 % goal** (Figure 1 C), judged on
  simulations not used to choose the demes: 98.5 % and 96.3 % for P36,
  90.5 % and 90.3 % for P34 (2025, 2026). Sampling every deme would
  reach 96 – 98 % at almost twice the cost.
- **The deme choice beats chance** (Figure 1 D): adding demes in random
  order until 90 % is reached takes a median of 42 mothers in P34 2025
  (the plan: 39) and 57 in P34 2026 (the plan: 45).
- **The choice is stable** (Figure 1 E): rerun on fresh simulations, most
  2025 demes are chosen every time; in 2026 a few near-equivalent plans
  exist, all reaching the goal, and the plan is fixed so the lab list
  does not change.

<a id="fig-1"></a>
![Figure 1](figures/step30_D_model_design_simulations.png)

**Figure 1.** Simulations behind the D1 plan. **A** — share of a deme's SRK alleles recovered against the number of mothers (15 seedlings each), one curve per deme size; dot = the number of mothers used; dotted line = 90 %. **B** — for each sampled deme, SRK alleles it holds (grey) and alleles its mothers' seedlings recover (red), mean and 95 % interval; labels give plants (pl.) and mothers (m.). **C** — share of the population's SRK alleles recovered by the plan (red) and by sampling every deme (blue), over independent simulations; bar = median. **D** — mothers needed to reach 90 % when demes are added in random order (grey) against the plan (red line). **E** — how often each deme is selected when the plan is rerun on 20 fresh simulations; red = in the plan.

---

## 5. Sampling

Running the design (§ 4) on the demes and mothers collected in 2025 and
2026 gives the plan below and the lab sheet that implements it.

### 5.1 The D1 plan

| Population | Season | Demes sampled (of usable) | Mothers | Seeds to germinate → seedlings | SRK alleles recovered (expected) |
|---|:-:|:-:|:-:|:-:|:-:|
| P36 (small) | 2025 | 2 of 2 | 4 | 100 → ~60 | 98.5 % |
| | 2026 | 2 of 2 | 3 | 75 → ~45 | 96.3 % |
| P34 (big) | 2025 | 9 of 14 | 39 | 975 → ~585 | 90.5 % |
| | 2026 | 4 of 20 | 45 | 1 110 → ~666 | 90.3 % |
| **Total** | | | **91** | **2 260 → ~1 356** | |

- **Sampled demes** (deme · mothers). P36: demes 1 and 2 in both
  seasons. P34 2025: 2 · 8, 3 · 6, 4 · 2, 5 · 9 (the adjusted deme,
  § 4.4), 8 · 2, 10 · 1, 12 · 4, 13 · 3, 14 · 4. P34 2026: 2 · 9, 3 · 5,
  9 · 17, 20 · 14.
- The sampled demes hold 76 % (2025) and 54 % (2026) of P34's fertile
  plants and span the population from north-west to south.
- **All 91 mothers have cleaned seed ready to germinate.** One P34 2026
  lot holds 10 seeds and is germinated in full.

<a id="fig-2"></a>
![Figure 2 — P36](figures/step30_D_model_design_demes_schematic_P36.png)

**Figure 2.** P36 (small) demes, 2025 (left) and 2026 (right). Coloured = sampled (label: deme · mothers); grey = not sampled (label: deme); dashed = a single-plant deme, left out. Black dots = slick spots, sized by fertile plants. Positions are in metres from an arbitrary origin, north up; the map carries no geographic coordinates. Both demes are sampled in both seasons. File: [`figures/step30_D_model_design_demes_schematic_P36.png`](figures/step30_D_model_design_demes_schematic_P36.png) ([PDF](figures/step30_D_model_design_demes_schematic_P36.pdf)).

<a id="fig-3"></a>
![Figure 3 — P34](figures/step30_D_model_design_demes_schematic_P34.png)

**Figure 3.** P34 (big) demes, 2025 (left) and 2026 (right); same conventions as Figure 2; * = the 2025 deme sampled with fewer mothers than its target (§ 4.4). File: [`figures/step30_D_model_design_demes_schematic_P34.png`](figures/step30_D_model_design_demes_schematic_P34.png) ([PDF](figures/step30_D_model_design_demes_schematic_P34.pdf)).

### 5.2 The lab sheet

**[`tables/step30_D_model_design_lab_sheet.tsv`](tables/step30_D_model_design_lab_sheet.tsv)**
— one row per mother (91 rows), sorted by population (small P36 first, then big P34), season, deme
and slick spot. It is keyed on **occurrenceID** (the wild mother), so it
works even where a germplasmID is still missing.

| Column | Meaning |
|---|---|
| `occurrenceID` | The wild mother plant — the key for the lab |
| `germplasmID` | Her cleaned seed lot |
| `BL`, `population`, `EO`, `locationCode` | Where she grows |
| `year` | Season the seed was collected (2025 or 2026) |
| `demeID`, `eventID_slick_spot` | Her deme and slick spot (Figures 2 and 3) |
| `pair`, `role` | D1, small (P36) or big (P34) |
| `status` | `germinate` — the seed lot is ready |
| `seeds_available` | Seeds in the lot |
| `seeds_to_germinate` | 25, or the whole lot if smaller |
| `seedlings_expected` | Seedlings to genotype (about 15) |
| `selected_because` | How she was drawn (one per slick spot, or among the deme's other mothers) |

**Per mother:** germinate `seeds_to_germinate` seeds, genotype the
seedlings (about 15), and record each seedling against the mother's
occurrenceID and germplasmID, so seedlings can be grouped by mother and
deme.

**Supporting files.**
- [`tables/step30_D_model_design_germplasm_by_deme.csv`](tables/step30_D_model_design_germplasm_by_deme.csv)
  — the same mothers listed by germplasmID, deme by deme (for pulling
  seed lots).
- [`tables/step30_D_model_design_demes.tsv`](tables/step30_D_model_design_demes.tsv)
  — every deme of P36 and P34: plants, mothers needed and held, whether
  sampled.
- [`tables/step30_D_model_design_summary.tsv`](tables/step30_D_model_design_summary.tsv)
  — the totals of § 5.1.

The plan is produced by `Scripts/step28_model_population_design.py`
(`./run_pipeline.sh design`); rerunning it gives the same lists.

---

## 6. Beyond D1, and what comes next

**Other populations.** With the mothers already collected, 27 of 33
populations (2025) and 32 of 37 (2026) hold enough mothers in every deme
for their seedlings to recover 90 % of their SRK alleles
([Figure 4](#fig-4)); the others need more mothers before their SRK
diversity can be measured this way.

<a id="fig-4"></a>
![Figure 4](figures/step30_A_sampling_adequacy.png)

**Figure 4.** Can each population's SRK diversity be measured with the mothers already collected? One row per population, grouped by Bottleneck Lineage; open circle = 2025, filled = 2026; red = not yet. **A** — expected share of the population's SRK alleles its collected mothers' seedlings would recover (dotted line = 90 %). **B** — share of fertile plants in demes with no collected mother. **C** — mothers collected (circle), mothers needed (|) and mothers with cleaned seed so far (×). Data: [`step30_A_sampling_adequacy_2025.tsv`](tables/step30_A_sampling_adequacy_2025.tsv), [`_2026`](tables/step30_A_sampling_adequacy_2026.tsv).

**Next: sampling for pollen compatibility.** The D1 genotypes will also
let us check the project's predictions of pollen compatibility — how
many plants each plant can mate with — at P36 and P34. Testing whether
plants with fewer compatible mates set less seed (mate limitation) is a
different question with different needs: a first simulation suggests it
takes many more mothers than D1, drawn across populations
([`figures/step23_power.png`](figures/step23_power.png)). That sampling
will be designed as a separate step, building on the SRK diversity
measured here.

**Assumptions to keep in mind.** The simulations treat each deme as an
independent draw from the species-wide allele frequencies. Real demes
have drifted further and neighbouring demes share alleles, so the real
populations probably hold fewer alleles than simulated — which makes the
90 % goal easier, not harder, to reach. Genotyped plants are treated as
independent; siblings of one mother are not yet corrected for.
