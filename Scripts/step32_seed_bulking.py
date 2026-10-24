"""step32_seed_bulking.py — reproductive status of every population, for
choosing source populations for seed bulking (talk + Seed_bulking.md).

Every measured plant enters the seed-by-size curve, but only populations
with seed data in **both seasons** (≥ `MIN_PLANTS` plants each in 2025 and
2026) are typed (Sven, 2026-10-24). Their type is the Step 13 call that
allows for the 2026 drop shared by all populations
(`priority_year_adjusted`):

- **Mate-limited** — sets less seed than its plants' size predicts
  (`mate_limited`; `candidate_mate_limited` drawn hollow);
- **Over-producing** — sets more seed than its size predicts, a possible
  sign of SI escape (`si_escape_candidate`);
- **Reproductively intact** — sets the seed its size predicts
  (`reproduction_intact`); these are the model populations for bulking;
- **Not yet clear** — switching between seasons or too few plants;
- **One season only** — populations with seed data in one season, shown
  in light grey and not typed.

The three type colours (magenta / gold / teal, grey for not yet clear)
are used on every figure and map of the talk and differ from the BL
colours. They are written to the population table so the R maps
(`step32_seed_bulking_maps.R`) read the same colours.

Figures (talk style, no BL colours) → figures_seed_bulking/ (local):
  seed_by_size_points        the same frame, points only (animation step 1)
  seed_by_size_plants        every plant, seeds vs crown size; band between
                             half and twice the seed its size predicts in
                             the type colours (magenta / teal / gold)
  population_seed_yield_frame  the same frame without data (animation step 1)
  population_seed_yield      seed yield for size per population and season,
                             coloured by type
  plant_size_examples        2026 field photos of a small, a medium and a
                             large plant on the measuring board (location
                             metadata stripped; copies in plant_photos/)
  genetic_example_P42_P34    SRK genotypes of a mate-limited (P42) and an
                             intact (P34) population: allele make-up,
                             mates available per plant, seed set for size
  bulking_filter_funnel      44 populations → seed measured → two seasons
                             (by type) → reproductively intact
  bulking_seed_available     cleaned seed held per candidate population and
                             season, with the number of mothers
Table: tables/step32_seed_bulking_populations.tsv (no coordinates)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from srk_units import load_event_units
from step13_size_seed_yield import load_occurrences

TABLES = Path("tables")
OUT = Path("figures_seed_bulking")
TYPE_OF = {"mate_limited": "Mate-limited",
           "si_escape_candidate": "Over-producing",
           "reproduction_intact": "Reproductively intact"}
TYPE_OF["candidate_mate_limited"] = "Mate-limited"   # drawn hollow
NOT_CLEAR = "Not yet clear"
ONE_SEASON = "One season only"
MIN_PLANTS = 3
ALPHA = 0.05        # season significant when Step 13's year-adjusted FDR q < ALPHA
TYPE_ORDER = ["Mate-limited", "Reproductively intact", "Over-producing", NOT_CLEAR, ONE_SEASON]
TYPE_COLORS = {"Mate-limited": "#C2185B",          # magenta
               "Over-producing": "#E6A100",        # gold
               "Reproductively intact": "#00897B",  # teal
               NOT_CLEAR: "#9E9E9E",               # grey
               ONE_SEASON: "#D9D9D9"}              # light grey
BAND = 2.0          # lines at twice and half the seed predicted by size
plt.rcParams.update({"font.size": 15, "axes.titlesize": 17, "axes.labelsize": 16,
                     "legend.fontsize": 13, "xtick.labelsize": 13, "ytick.labelsize": 13})


def population_table() -> pd.DataFrame:
    c = pd.read_csv(TABLES / "step30i_conservation_priority.tsv", sep="\t")
    c["two_seasons"] = ((c["n_plants_2025"].fillna(0) >= MIN_PLANTS)
                        & (c["n_plants_2026"].fillna(0) >= MIN_PLANTS))
    c["type"] = np.where(c["two_seasons"],
                         c["priority_year_adjusted"].map(TYPE_OF).fillna(NOT_CLEAR), ONE_SEASON)
    c["candidate"] = c["two_seasons"] & (c["priority_year_adjusted"] == "candidate_mate_limited")
    c["color"] = c["type"].map(TYPE_COLORS)
    c["model_for_bulking"] = c["type"] == "Reproductively intact"
    cl = pd.read_csv(TABLES / "step30g_populations_classified.tsv", sep="\t")
    c = c.merge(cl[["populationID", "locationCodes"]], on="populationID", how="left")
    plants = []
    for yr in (2025, 2026):
        d = pd.read_csv(TABLES / f"step30_A_sampling_adequacy_per_deme_{yr}.tsv", sep="\t")
        g = d.groupby("populationID").agg(**{f"fertile_plants_{yr}": ("deme_N_fertile", "sum"),
                                              f"mothers_collected_{yr}": ("n_occurrences", "sum")})
        plants.append(g)
    c = c.merge(pd.concat(plants, axis=1).reset_index(), on="populationID", how="left")
    s = pd.read_csv(TABLES / "step30i_size_seed_population_strata.tsv", sep="\t")
    q = s.pivot(index="populationID", columns="year", values="q_year_adjusted")
    q.columns = [f"q_{y}" for y in q.columns]
    fa = s.pivot(index="populationID", columns="year", values="fold_year_adjusted")
    fa.columns = [f"fold_adj_{y}" for y in fa.columns]
    c = c.merge(q.reset_index(), on="populationID", how="left")
    c = c.merge(fa.reset_index(), on="populationID", how="left")
    cols = ["populationID", "BL", "locationCodes", "two_seasons", "type", "candidate", "color",
            "model_for_bulking", "priority_year_adjusted", "robust_priority", "fold_2025", "fold_2026", "n_plants_2025", "n_plants_2026",
            "fold_adj_2025", "fold_adj_2026", "q_2025", "q_2026", "fertile_plants_2025", "fertile_plants_2026",
            "mothers_collected_2025", "mothers_collected_2026"]
    t = c[cols].copy()
    t["type"] = pd.Categorical(t["type"], TYPE_ORDER, ordered=True)
    return t.sort_values(["type", "BL", "populationID"]).reset_index(drop=True)


def plants_with_type(pops: pd.DataFrame) -> pd.DataFrame:
    occ = load_occurrences()
    u = pd.concat([load_event_units(y)[["eventID", "populationID"]].assign(year=y)
                   for y in (2025, 2026)])
    m = occ.merge(u, on=["eventID", "year"]).merge(pops[["populationID", "type"]], on="populationID")
    return m


def seed_by_size(m: pd.DataFrame, name: str, bands: bool = True) -> None:
    """Every plant in one colour; the band between half and twice the seed
    a plant's size predicts shaded teal (as expected), the upper line gold
    (more seed) and the lower line magenta (less seed) — the colour codes
    of the reproductive types, leading into population_seed_yield.
    bands=False draws the same frame with the points only (band, lines
    and labels invisible), so the two images overlay exactly in a slide
    animation."""
    show = 1.0 if bands else 0.0
    b1, b0 = np.polyfit(m["logc"], m["logy"], 1)          # pooled size → seed curve
    fig, ax = plt.subplots(figsize=(11, 7.5))
    xs = np.linspace(m["logc"].min(), m["logc"].max(), 100)
    lo = 10 ** (b0 + b1 * xs + np.log10(1 / BAND))
    hi = 10 ** (b0 + b1 * xs + np.log10(BAND))
    ax.fill_between(10 ** xs, lo, hi, color=TYPE_COLORS["Reproductively intact"],
                    alpha=0.16 * show, lw=0, zorder=1)
    ax.scatter(10 ** m["logc"], 10 ** m["logy"], s=22, color="#4D4D4D", alpha=0.5,
               linewidths=0, zorder=2)
    ax.plot(10 ** xs, hi, color=TYPE_COLORS["Over-producing"], lw=3, zorder=3, alpha=show)
    ax.plot(10 ** xs, lo, color=TYPE_COLORS["Mate-limited"], lw=3, zorder=3, alpha=show)
    mid = 10 ** (b0 + b1 * xs[-1])
    for y, lab, tp in ((hi[-1], "↑ over-producing\n(more seed, ×2)", "Over-producing"),
                       (mid, "as expected", "Reproductively intact"),
                       (lo[-1], "↓ mate-limited\n(less seed, ÷2)", "Mate-limited")):
        ax.text(10 ** xs[-1] * 1.04, y, lab, va="center", ha="left", fontsize=14,
                weight="bold", color=TYPE_COLORS[tp], alpha=show)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(10 ** m["logc"].min() * 0.9, 10 ** m["logc"].max() * 2.2)
    ax.set_ylim(10 ** m["logy"].min() * 0.8, 10 ** m["logy"].max() * 1.25)
    ax.set_xlabel("Plant size (crown diameter, cm)")
    ax.set_ylabel("Seeds per plant")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.01, 0.99, f"{len(m):,} plants · 2025 + 2026", transform=ax.transAxes,
            va="top", fontsize=13, color="#444")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def population_seed_yield(pops: pd.DataFrame, data: bool = True,
                          name: str = "population_seed_yield") -> None:
    """Fold of expectation per season for the two-season populations,
    ordered within type; one-season populations as one light-grey row.
    data=False draws the same frame with the data invisible (axes, labels,
    lines and legend kept), for a slide animation."""
    show = 1.0 if data else 0.0
    two = pops[pops["two_seasons"]].copy()
    two["keyf"] = np.exp((np.log(two["fold_adj_2025"]) + np.log(two["fold_adj_2026"])) / 2)
    two = two.sort_values(["type", "keyf"], ascending=[False, True])
    one = pops[~pops["two_seasons"]]
    n = len(two)
    fig, ax = plt.subplots(figsize=(11, 0.55 * (n + 1) + 2.2))
    for y, (_, r) in enumerate(two.iterrows(), start=1):
        c = TYPE_COLORS[r["type"]]
        dy = {2025: 0.13, 2026: -0.13}          # seasons slightly apart so both stay visible
        ax.plot([r["fold_adj_2025"], r["fold_adj_2026"]], [y + dy[2025], y + dy[2026]],
                color=c, lw=2, alpha=0.7 * show)
        for yr, mk in ((2025, "o"), (2026, "s")):
            sig = pd.notna(r[f"q_{yr}"]) and r[f"q_{yr}"] < ALPHA
            ax.scatter(r[f"fold_adj_{yr}"], y + dy[yr], s=130, marker=mk,
                       facecolors=c if sig else "white", edgecolors=c, linewidths=2.2, zorder=3,
                       alpha=show)
    f1 = pd.concat([one["fold_adj_2025"], one["fold_adj_2026"]]).dropna()
    ax.scatter(f1, np.zeros(len(f1)), s=45, color=TYPE_COLORS[ONE_SEASON], zorder=2, alpha=show)
    labels = [ONE_SEASON + f" ({len(one)})"] + [
        f"P{r['populationID']} ({r['locationCodes'].split(',')[0].strip()})" for _, r in two.iterrows()]
    ax.set_yticks(range(n + 1)); ax.set_yticklabels(labels, fontsize=13)
    ax.get_yticklabels()[0].set_color("#888888")
    ax.axvline(1, color="black", lw=1.2)
    for f in (BAND, 1 / BAND):
        ax.axvline(f, color="black", lw=1, ls="--")
    ax.set_xscale("log")
    ax.set_xticks([0.02, 0.1, 0.25, 0.5, 1, 2, 4])
    ax.set_xticklabels(["0.02", "0.1", "0.25", "0.5", "1", "2", "4"])
    ax.set_xlabel("Seeds per plant ÷ seeds expected for its size\n"
                  "(allowing for the 2026 drop shared by all populations)")
    for xv, yv, lab, tp, ha in ((1 / BAND, 1.01, "← mate-limited", "Mate-limited", "right"),
                                (1, 1.045, "as expected", "Reproductively intact", "center"),
                                (BAND, 1.01, "over-producing →", "Over-producing", "left")):
        ax.text(xv, yv, lab, transform=ax.get_xaxis_transform(), ha=ha, va="bottom",
                fontsize=13, weight="bold", color=TYPE_COLORS[tp])
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=TYPE_COLORS[t], ms=11, label=t)
         for t in ("Mate-limited", "Reproductively intact", "Over-producing", NOT_CLEAR)]
    h += [Line2D([], [], marker="o", ls="", mfc="#444", mec="#444", ms=11, label="2025"),
          Line2D([], [], marker="s", ls="", mfc="#444", mec="#444", ms=11, label="2026"),
          Line2D([], [], marker="o", ls="", mfc="#444", mec="#444", ms=11,
                 label="filled = differs from expectation\n(q < 0.05)"),
          Line2D([], [], marker="o", ls="", mfc="white", mec="#444", mew=2, ms=11,
                 label="open = not significant")]
    ax.legend(handles=h, loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=12, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def filter_funnel(pops: pd.DataFrame) -> None:
    """Horizontal funnel: all populations → seed measured → two seasons
    (split by type) → reproductively intact → genetic filters (to come)."""
    n_all = len(pd.read_csv(TABLES / "step30g_populations_classified.tsv", sep="\t"))
    two = pops[pops["two_seasons"]]
    counts = two["type"].value_counts()
    stages = [("All populations", n_all, None),
              ("Seed set measured", len(pops), None),
              ("Measured in both seasons", len(two), counts),
              ("Reproductively intact", int(counts.get("Reproductively intact", 0)), None)]
    fig, ax = plt.subplots(figsize=(12, 5.8))
    y = np.arange(len(stages))[::-1]
    for yi, (lab, n, split) in zip(y, stages):
        if split is None:
            col = TYPE_COLORS["Reproductively intact"] if lab.startswith("Reproductively") else "#B0B0B0"
            ax.barh(yi, n, color=col, height=0.62)
        else:
            left = 0
            for tp in ("Mate-limited", "Reproductively intact", "Over-producing", NOT_CLEAR):
                k = int(split.get(tp, 0))
                ax.barh(yi, k, left=left, color=TYPE_COLORS[tp], height=0.62)
                if k:
                    ax.text(left + k / 2, yi, str(k), ha="center", va="center", color="white",
                            fontsize=14, weight="bold")
                left += k
        ax.text(n + 0.6, yi, f"{n}", va="center", fontsize=17, weight="bold")
    ax.set_yticks(y); ax.set_yticklabels([s[0] for s in stages], fontsize=16)
    ax.text(0.6, -0.95, "Next: genetic filters — SRK diversity, self-compatibility, one source per lineage",
            fontsize=14, style="italic", color="#444")
    ax.set_xlim(0, n_all * 1.12); ax.set_ylim(-1.3, len(stages) - 0.5)
    ax.set_xlabel("Populations")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=TYPE_COLORS[tp], label=tp) for tp in
                       ("Mate-limited", "Reproductively intact", "Over-producing", NOT_CLEAR)],
              loc="lower right", bbox_to_anchor=(1.0, 0.18), fontsize=13, frameon=False,
              title="Measured in both seasons", title_fontsize=13)
    ax.spines[["top", "right", "left"]].set_visible(False); ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"bulking_filter_funnel.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def seed_available(pops: pd.DataFrame) -> pd.DataFrame:
    """Cleaned seed held per candidate population and season (LEPA DB, as in
    the sampling design), with mothers collected and lots still to clean."""
    from step28_model_population_design import load_mothers
    m = load_mothers()
    cand = pops.loc[pops["model_for_bulking"], "populationID"].tolist()
    rows = []
    for yr in (2025, 2026):
        x = m.merge(load_event_units(yr)[["eventID", "populationID"]], on="eventID")
        for p in cand:
            g = x[x["populationID"] == p]
            c = g[g["cleaned"] == 1]
            rows.append(dict(populationID=p, year=yr, mothers=g["occurrenceID"].nunique(),
                             lots_cleaned=len(c), lots_to_clean=len(g) - len(c),
                             seeds_cleaned=int(round(c["seeds"].sum()))))
    s = pd.DataFrame(rows)
    s.to_csv(TABLES / "step32_seed_available.tsv", sep="\t", index=False)
    lab = {r["populationID"]: f"P{r['populationID']}\n({r['locationCodes'].split(',')[0].strip()})"
           for _, r in pops[pops["model_for_bulking"]].iterrows()}
    fig, ax = plt.subplots(figsize=(11, 6.5))
    xpos = np.arange(len(cand)); w = 0.38
    teal = TYPE_COLORS["Reproductively intact"]
    for k, yr in enumerate((2025, 2026)):
        g = s[s["year"] == yr].set_index("populationID").loc[cand]
        bars = ax.bar(xpos + (k - 0.5) * w, g["seeds_cleaned"] / 1000, width=w,
                      color=teal if yr == 2026 else "white", edgecolor=teal, linewidth=2.5,
                      label=str(yr))
        for b, (_, r) in zip(bars, g.iterrows()):
            note = f"{r['mothers']} mothers"
            if r["lots_to_clean"]:
                note += f"\n({r['lots_to_clean']} lots\nto clean)"
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 2, note, ha="center",
                    va="bottom", fontsize=11)
    ax.set_xticks(xpos); ax.set_xticklabels([lab[p] for p in cand], fontsize=15)
    ax.set_ylabel("Cleaned seed held (thousands)")
    ax.set_ylim(0, s["seeds_cleaned"].max() / 1000 * 1.3)
    ax.legend(fontsize=13, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"bulking_seed_available.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    return s


SIZE_EXAMPLES = [("Small", 2515), ("Medium", 2592), ("Large", 2748)]   # occurrenceIDs, 2026


def plant_size_panel() -> None:
    """Three 2026 measuring-board photos (small / medium / large plant, each
    near the seed its size predicts), saved without any metadata."""
    import sqlite3
    from PIL import Image, ImageOps
    from project_paths import LEPA_DB
    con = sqlite3.connect(f"file:{LEPA_DB}?mode=ro", uri=True)
    ids = ",".join(str(o) for _, o in SIZE_EXAMPLES)
    d = pd.read_sql(f"""SELECT p.occurrenceID, p.occurrenceCrownSize AS crown, m.identifier,
                               t.seedQuantityTotal AS seeds
                        FROM Phenotyping p JOIN Multimedia m ON m.multimediaID = p.multimediaID
                        JOIN vOccurrenceTraits t ON t.occurrenceID = p.occurrenceID
                        WHERE p.occurrenceID IN ({ids})""", con).set_index("occurrenceID")
    con.close()
    photos = Path(LEPA_DB).parent / "Multimedia_main"
    out = OUT / "plant_photos"
    out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.6))
    for ax, (lab, oc) in zip(axes, SIZE_EXAMPLES):
        r = d.loc[oc]
        im = ImageOps.exif_transpose(Image.open(photos / r["identifier"])).convert("RGB")
        im.thumbnail((2400, 2400))
        clean = Image.new("RGB", im.size)
        clean.putdata(list(im.getdata()))          # pixels only: no GPS or other metadata
        clean.save(out / f"size_{lab.lower()}_OC{oc}_crown{r['crown']:g}cm_{r['seeds']:.0f}seeds.jpg",
                   quality=90)
        ax.imshow(clean); ax.axis("off")
        ax.set_title(f"{lab} plant", fontsize=20, weight="bold", pad=8)
        ax.text(0.5, -0.06, f"crown {r['crown']:g} cm  ·  {r['seeds']:,.0f} seeds",
                transform=ax.transAxes, ha="center", va="top", fontsize=17)
    fig.tight_layout(w_pad=1.5)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"plant_size_examples.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def genetic_example(pops: pd.DataFrame) -> pd.DataFrame:
    """P42 (mate-limited) against P34 (reproductively intact) with their SRK
    genotypes, the plants and compatibility statistic of Step 25: each
    plant against fathers drawn from its population's observed allele
    frequencies, under sporophytic Class I / II SI with the empirical
    zygosity. Alleles are labelled A, B, … by their frequency over the two
    populations (codes in the TSV)."""
    import step25_partC_clean_overlap as s25
    from step24_EO_scale_validation import compute_local_f, load_allele_to_fg, load_individuals
    from srk_si_model import (build_class_i_mask, load_class_map, load_zygosity_dist,
                              p_compat_sporophytic_empirical)
    rng = np.random.default_rng(2038)
    prior = pd.read_csv(s25.DEFAULT_PRIOR_TSV, sep="\t", encoding="utf-8-sig")
    fg = prior["Fg"].astype(str).tolist()
    ind = load_individuals(s25.DEFAULT_INDIV_TSV, load_allele_to_fg(s25.DEFAULT_ALLELE_FG_TSV), fg)
    pred_pc = pd.read_csv(s25.DEFAULT_PRED_PCOMPAT, sep="\t", encoding="utf-8-sig")
    pred_dv = pd.read_csv(s25.DEFAULT_PRED_DIVERSITY, sep="\t", encoding="utf-8-sig")
    ind, _, _, _ = s25.add_population_sites(ind, pred_pc, pred_dv)
    class_i = build_class_i_mask(fg, load_class_map())
    zyg = load_zygosity_dist()
    sites = {42: "EO70", 34: "P34"}
    res, freqs = {}, {}
    for pop, site in sites.items():
        sub = ind[(ind["EO"] == site) & (ind["n_functional_copies"] > 0)]
        g = np.stack(sub["genotype"].values)
        f = compute_local_f(np.stack(sub["fg_copies"].values))
        comp = p_compat_sporophytic_empirical(g, f, class_i, zyg, n_fathers=2000, rng=rng)
        res[pop] = dict(n=len(sub), compat=comp, alleles=int((f > 0).sum()))
        freqs[pop] = f
    order = np.argsort(-(freqs[42] + freqs[34]))
    order = [i for i in order if freqs[42][i] + freqs[34][i] > 0]
    letter = {i: chr(65 + k) for k, i in enumerate(order)}
    rows = [dict(allele=letter[i], srk_group=fg[i], **{f"freq_P{p}": float(freqs[p][i]) for p in sites})
            for i in order]
    tab = pd.DataFrame(rows)
    tab.to_csv(TABLES / "step32_genetic_example_P42_P34.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(18, 7))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1, 0.85], wspace=0.55)
    col = {42: TYPE_COLORS["Mate-limited"], 34: TYPE_COLORS["Reproductively intact"]}
    name = {42: "P42 (EO70)\nmate-limited", 34: "P34 (EO18)\nreproductively intact"}
    # A — allele make-up
    ax = fig.add_subplot(gs[0])
    shades = plt.get_cmap("Greys")(np.linspace(0.85, 0.25, len(order)))
    for yi, pop in zip((1, 0), (42, 34)):
        left = 0.0
        for k, i in enumerate(order):
            w = freqs[pop][i]
            if w <= 0:
                continue
            ax.barh(yi, w, left=left, color=shades[k], edgecolor="white", height=0.6)
            if w >= 0.06:
                ax.text(left + w / 2, yi, letter[i], ha="center", va="center", fontsize=14,
                        weight="bold", color="white" if k < len(order) / 2 else "black")
            left += w
        top2 = np.sort(freqs[pop])[::-1][:2].sum()
        ax.text(0, yi + 0.34, f"{res[pop]['alleles']} alleles · top two hold {100 * top2:.0f} %",
                va="bottom", ha="left", fontsize=13, color=col[pop], weight="bold")
    ax.set_ylim(-0.45, 1.6)
    ax.set_yticks([1, 0]); ax.set_yticklabels([name[42], name[34]], fontsize=14)
    for lab, pop in zip(ax.get_yticklabels(), (42, 34)):
        lab.set_color(col[pop])
    from matplotlib.ticker import PercentFormatter
    ax.set_xlim(0, 1); ax.set_xlabel("Share of SRK allele copies")
    ax.xaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.set_title("A. SRK allele make-up", loc="left", weight="bold")
    # B — mates available per plant
    ax = fig.add_subplot(gs[1])
    for xi, pop in zip((0, 1), (42, 34)):
        c = res[pop]["compat"]
        ax.scatter(xi + rng.uniform(-0.18, 0.18, len(c)), c, s=30, color=col[pop], alpha=0.6,
                   linewidths=0)
        ax.plot([xi - 0.3, xi + 0.3], [c.mean()] * 2, color="black", lw=3)
        ax.text(xi + 0.33, c.mean(), f"{100 * c.mean():.0f} %", va="center", fontsize=15,
                weight="bold")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["P42", "P34"], fontsize=15)
    for lab, pop in zip(ax.get_xticklabels(), (42, 34)):
        lab.set_color(col[pop])
    ax.set_ylim(0, 1.02); ax.set_xlim(-0.6, 1.8)
    from matplotlib.ticker import PercentFormatter
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.set_ylabel("Share of the population\na plant can mate with")
    ax.set_title("B. Mates available per plant", loc="left", weight="bold")
    # C — seed set for size
    ax = fig.add_subplot(gs[2])
    pt = pops.set_index("populationID")
    for xi, pop in zip((0, 1), (42, 34)):
        r = pt.loc[pop]
        for yr, mk, dx in ((2025, "o", -0.12), (2026, "s", 0.12)):
            sig = pd.notna(r[f"q_{yr}"]) and r[f"q_{yr}"] < ALPHA
            ax.scatter(xi + dx, r[f"fold_adj_{yr}"], s=170, marker=mk, edgecolors=col[pop],
                       facecolors=col[pop] if sig else "white", linewidths=2.4, zorder=3)
    ax.axhline(1, color="black", lw=1.2)
    for fb in (BAND, 1 / BAND):
        ax.axhline(fb, color="black", lw=1, ls="--")
    ax.set_yscale("log"); ax.set_ylim(0.35, 3.6)
    ax.set_yticks([0.5, 1, 2]); ax.set_yticklabels(["0.5", "1", "2"])
    ax.yaxis.set_minor_locator(plt.NullLocator())
    ax.set_xticks([0, 1]); ax.set_xticklabels(["P42", "P34"], fontsize=15)
    for lab, pop in zip(ax.get_xticklabels(), (42, 34)):
        lab.set_color(col[pop])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylabel("Seeds per plant ÷ seeds\nexpected for its size")
    for yv, lab, tp in ((BAND, "over-producing", "Over-producing"),
                        (1, "as expected", "Reproductively intact"),
                        (1 / BAND, "mate-limited", "Mate-limited")):
        ax.text(1.02, yv, lab, transform=ax.get_yaxis_transform(), va="center", ha="left",
                fontsize=13, weight="bold", color=TYPE_COLORS[tp])
    ax.set_title("C. Seed set for plant size", loc="left", weight="bold")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color="#444", ms=10, label="2025"),
                       Line2D([], [], marker="s", ls="", color="#444", ms=10, label="2026"),
                       Line2D([], [], marker="o", ls="", mfc="white", mec="#444", mew=2, ms=10,
                              label="open = not significant")],
              fontsize=11, loc="upper left", bbox_to_anchor=(0, 1.0), frameon=False)
    for a in fig.axes:
        a.spines[["top", "right"]].set_visible(False)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"genetic_example_P42_P34.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[step32] genetic example:", {f"P{p}": dict(n=v["n"], alleles=v["alleles"],
          mean_compat=round(float(v["compat"].mean()), 3)) for p, v in res.items()})
    return tab


def main() -> None:
    OUT.mkdir(exist_ok=True)
    pops = population_table()
    pops.to_csv(TABLES / "step32_seed_bulking_populations.tsv", sep="\t", index=False)
    print(pops.groupby("type", observed=True)["populationID"].apply(list).to_string())
    m = plants_with_type(pops)
    plant_size_panel()
    seed_by_size(m, "seed_by_size_plants")
    seed_by_size(m, "seed_by_size_points", bands=False)       # animation step 1
    population_seed_yield(pops)
    population_seed_yield(pops, data=False, name="population_seed_yield_frame")   # animation step 1
    filter_funnel(pops)
    genetic_example(pops)
    print(seed_available(pops).to_string(index=False))
    print(f"[step32] wrote tables/step32_seed_bulking_populations.tsv and figures in {OUT}/")


if __name__ == "__main__":
    main()
