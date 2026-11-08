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
                             mates available per plant, seed set for size;
                             talk build-up _step1 (headers: S-types, plants)
                             and _step2 (+ genotypes, 3 commonest sets)
  bulking_filter_funnel      44 populations → seed measured → two seasons
                             (by type) → reproductively intact
  bulking_seed_available     cleaned seed held per candidate population and
                             season — LOCAL ONLY (gitignored): seed counts are
                             not shared
Tables: tables/step32_seed_bulking_populations.tsv and step32_materials_by_season.tsv
(locations, EOs, slick spots and plants per season; no coordinates)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from srk_units import load_event_units
from step13_size_seed_yield import load_occurrences

with plt.rc_context():                  # the cartoon modules set their own font size
    from type_icons import icon as type_icon


def add_type_icons(fig, texts, side: str, size_in: float = 0.95) -> None:
    """Badge icon of each reproductive type next to its label — `texts` maps
    type → matplotlib Text; side "above" (centred over the label) or "right"
    (after it). Icons from type_icons.py."""
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    w, h = size_in / fig.get_figwidth(), size_in / fig.get_figheight()
    for tp, t in texts.items():
        bb = t.get_window_extent()
        (x0, y0), (x1, y1) = inv.transform((bb.x0, bb.y0)), inv.transform((bb.x1, bb.y1))
        if side == "above":
            rect = [(x0 + x1) / 2 - w / 2, y1 + 0.004, w, h]
        else:
            rect = [x1 + 0.006, (y0 + y1) / 2 - h / 2, w, h]
        type_icon(fig.add_axes(rect), tp)

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


def season_significant(r, yr: int) -> bool:
    """Filled marker: for a reproductively intact population, seed set is
    significantly WITHIN half – twice what size predicts (Step 13's
    equivalence test); for any other type, it significantly DIFFERS from
    expectation."""
    q = r[f"qeq_{yr}"] if r["type"] == "Reproductively intact" else r[f"q_{yr}"]
    return bool(pd.notna(q) and q < ALPHA)
TYPE_ORDER = ["Mate-limited", "Reproductively intact", "Over-producing", NOT_CLEAR, ONE_SEASON]
TYPE_COLORS = {"Mate-limited": "#E6A100",          # gold
               "Over-producing": "#C2185B",        # magenta (close to the red "avoid" sign)
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
    qe = s.pivot(index="populationID", columns="year", values="q_equiv_year_adjusted")
    qe.columns = [f"qeq_{y}" for y in qe.columns]
    c = c.merge(q.reset_index(), on="populationID", how="left")
    c = c.merge(qe.reset_index(), on="populationID", how="left")
    c = c.merge(fa.reset_index(), on="populationID", how="left")
    cols = ["populationID", "BL", "locationCodes", "two_seasons", "type", "candidate", "color",
            "model_for_bulking", "priority_year_adjusted", "robust_priority", "fold_2025", "fold_2026", "n_plants_2025", "n_plants_2026",
            "fold_adj_2025", "fold_adj_2026", "q_2025", "q_2026", "qeq_2025", "qeq_2026",
            "fertile_plants_2025", "fertile_plants_2026",
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
    a plant's size predicts shaded teal (as expected), the fitted size →
    seed line teal through its middle, the upper line gold
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
    ax.plot(10 ** xs, 10 ** (b0 + b1 * xs), color=TYPE_COLORS["Reproductively intact"], lw=3,
            zorder=3, alpha=show)                                # the size → seed relationship
    mid = 10 ** (b0 + b1 * xs[-1])
    labels = {}
    for y, lab, tp in ((hi[-1], "↑ over-producing\n(more seed, ×2)", "Over-producing"),
                       (mid, "as expected", "Reproductively intact"),
                       (lo[-1], "↓ mate-limited\n(less seed, ÷2)", "Mate-limited")):
        labels[tp] = ax.text(10 ** xs[-1] * 1.04, y, lab, va="center", ha="left", fontsize=14,
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
    if bands:                                   # icons in the empty space each line bounds
        fig.canvas.draw()
        inv = fig.transFigure.inverted()
        size_in = 1.35
        w, h = size_in / fig.get_figwidth(), size_in / fig.get_figheight()
        for tp, (xd, yd) in (("Over-producing", (1.9, 2500)),     # above the gold line
                             ("Mate-limited", (40, 12))):         # below the magenta line
            cx, cy = inv.transform(ax.transData.transform((xd, yd)))
            type_icon(fig.add_axes([cx - w / 2, cy - h / 2, w, h]), tp)
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
            sig = season_significant(r, yr)
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
    ax.set_xlabel("Seeds per plant ÷ seeds expected for its size")
    labels = {}
    for xv, yv, lab, tp, ha in ((1 / BAND, 1.01, "← mate-limited", "Mate-limited", "right"),
                                (1, 1.045, "as expected", "Reproductively intact", "center"),
                                (BAND, 1.01, "over-producing →", "Over-producing", "left")):
        labels[tp] = ax.text(xv, yv, lab, transform=ax.get_xaxis_transform(), ha=ha,
                             va="bottom", fontsize=13, weight="bold", color=TYPE_COLORS[tp])
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=TYPE_COLORS[t], ms=11, label=t)
         for t in ("Mate-limited", "Reproductively intact", "Over-producing", NOT_CLEAR)]
    h += [Line2D([], [], marker="o", ls="", mfc="#444", mec="#444", ms=11, label="2025"),
          Line2D([], [], marker="s", ls="", mfc="#444", mec="#444", ms=11, label="2026"),
          Line2D([], [], marker="o", ls="", mfc="#444", mec="#444", ms=11,
                 label="filled = significant (q < 0.05): differs from\n"
                       "expectation, or (teal) lies within half – twice"),
          Line2D([], [], marker="o", ls="", mfc="white", mec="#444", mew=2, ms=9,
                 label="open = not significant"),
          Line2D([], [], ls="", label="x axis allows for the 2026 drop\n"
                                      "shared by all populations")]
    for hh in h:
        hh.set_markersize(9)
    leg = ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=3,
                    fontsize=10.5, frameon=True, edgecolor="#999999", fancybox=False,
                    title="Legend", title_fontproperties={"weight": "bold", "size": 11.5},
                    columnspacing=1.8, handletextpad=0.4, borderpad=0.7)
    leg.get_frame().set_linewidth(1.0)                     # boxed legend under the plot
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    add_type_icons(fig, labels, "above", size_in=1.25)
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
    ax.text(0.6, -0.95, "Next: genetic filters — S-allele diversity, self-compatibility, one source per lineage",
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
    s.to_csv(TABLES / "step32_seed_available.tsv", sep="\t", index=False)   # local only (gitignored)
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


def plant_size_panel(m: pd.DataFrame) -> None:
    """Three 2026 measuring-board photos (small / medium / large plant, each
    near the seed its size predicts), saved without any metadata; each is
    labelled with its crown size and the seed its size predicts on the
    Figure 9 curve (pooled size → seed fit over all plants). The crown shown
    is the one written on the field board in the photo (Phenotyping remarks
    "field board h=… w=…"), so label and board agree; the database value
    (measured afterwards from the photo's rulers) is used when no board
    value is recorded."""
    b1, b0 = np.polyfit(m["logc"], m["logy"], 1)
    import sqlite3
    from PIL import Image, ImageOps
    from project_paths import LEPA_DB
    con = sqlite3.connect(f"file:{LEPA_DB}?mode=ro", uri=True)
    ids = ",".join(str(o) for _, o in SIZE_EXAMPLES)
    d = pd.read_sql(f"""SELECT p.occurrenceID, p.occurrenceCrownSize AS crown, p.remarks, m.identifier,
                               t.seedQuantityTotal AS seeds
                        FROM Phenotyping p JOIN Multimedia m ON m.multimediaID = p.multimediaID
                        JOIN vOccurrenceTraits t ON t.occurrenceID = p.occurrenceID
                        WHERE p.occurrenceID IN ({ids})""", con).set_index("occurrenceID")
    con.close()
    board = d["remarks"].str.extract(r"board h=[\d.]+ w=([\d.]+)", expand=False).astype(float)
    d["crown"] = board.fillna(d["crown"])               # crown as written on the board
    photos = Path(LEPA_DB).parent / "Multimedia_main"
    out = OUT / "plant_photos"
    out.mkdir(parents=True, exist_ok=True)
    cleaned = {}
    for lab, oc in SIZE_EXAMPLES:
        r = d.loc[oc]
        im = ImageOps.exif_transpose(Image.open(photos / r["identifier"])).convert("RGB")
        im.thumbnail((2400, 2400))
        clean = Image.new("RGB", im.size)
        clean.putdata(list(im.getdata()))          # pixels only: no GPS or other metadata
        clean.save(out / f"size_{lab.lower()}_OC{oc}_crown{r['crown']:g}cm.jpg", quality=90)
        cleaned[oc] = clean

    def size_row(fig, axes):
        for ax, (lab, oc) in zip(axes, SIZE_EXAMPLES):
            crown = d.loc[oc, "crown"]
            ax.imshow(cleaned[oc]); ax.axis("off")
            ax.set_title(f"{lab} plant", fontsize=20, weight="bold", pad=8)
            exp_seed = 10 ** (b0 + b1 * np.log10(crown))
            ax.text(0.5, -0.06, f"crown {crown:g} cm",
                    transform=ax.transAxes, ha="center", va="top", fontsize=17)
            ax.text(0.5, -0.16, f"expected ≈ {round(exp_seed, -1):,.0f} seeds",
                    transform=ax.transAxes, ha="center", va="top", fontsize=17, weight="bold",
                    color=TYPE_COLORS["Reproductively intact"])

    # Figure 8: plant size (top) and seed yield (bottom)
    fig = plt.figure(figsize=(18, 14.4))
    fig.text(0.02, 0.985, "Plant size (field)", fontsize=22, weight="bold", va="top")
    size_row(fig, [fig.add_axes([0.0 + k * 0.337, 0.615, 0.326, 0.32]) for k in range(3)])
    fig.text(0.02, 0.52, "Seed yield (lab)", fontsize=22, weight="bold", va="top")
    seed_cleaning_panel(fig, out, bottom=0.0, height=0.44)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"plant_size_examples.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # talk: one 16:9 slide per row
    title = dict(fontsize=26, weight="bold", va="top")
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.03, 0.95, "Plant phenotyping: plant size (field)", **title)
    size_row(fig, [fig.add_axes([0.005 + k * 0.333, 0.32, 0.323, 0.52]) for k in range(3)])
    fig.text(0.5, 0.12, "Each plant is photographed against a board with rulers;"
             " its crown predicts the seed it should set", ha="center", fontsize=17, color="#444")
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"plant_size_examples_slide1.{ext}", dpi=200)
    plt.close(fig)
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.03, 0.95, "Plant phenotyping: seed yield (lab)", **title)
    seed_cleaning_panel(fig, out, bottom=0.08, height=0.62)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"plant_size_examples_slide2.{ext}", dpi=200)
    plt.close(fig)


# Seed-count model of the LEPA DB (germplasmQuantityEstimate): lot weight (g)
# = SEED_W0 + SEED_W1 x seeds, 95 % prediction interval +/- SEED_PI g
SEED_W0, SEED_W1, SEED_PI = 0.0003681555, 0.0004252451, 0.01926434


def seed_cleaning_panel(fig, out: Path, bottom: float, height: float) -> None:
    """Seed-yield row of Figure 8 (and talk slide 2): the cleaning bench for
    context, then the close-up of a cleaned lot (seed separated from chaff)
    with the weighing photo inset at its top right and, under it, the
    weight → seed-count model as a legend. Both photos at the same height,
    centred; they are metadata-free copies (plant_photos/
    seed_cleaning_lab.jpg, seed_weighing_lab.jpg), cropped here."""
    from PIL import Image
    src = Image.open(out / "seed_cleaning_lab.jpg")
    bench = src.crop((0, 375, 1440, 1650))
    close = src.crop((705, 960, 1800, 1665))
    fw, fh = fig.get_size_inches()
    h_in, gap = height * fh, 0.25
    w1, w2 = h_in * bench.width / bench.height, h_in * close.width / close.height
    x = (fw - w1 - w2 - gap) / 2
    title = dict(fontsize=20, weight="bold", pad=8)
    ax1 = fig.add_axes([x / fw, bottom, w1 / fw, height])
    ax1.imshow(bench); ax1.axis("off"); ax1.set_title("Cleaning by hand", **title)
    ax = fig.add_axes([(x + w1 + gap) / fw, bottom, w2 / fw, height])
    ax.imshow(close); ax.axis("off"); ax.set_title("Cleaned seed, weighed and counted", **title)
    ax.text(0.03, 0.64, "chaff", transform=ax.transAxes, fontsize=17, weight="bold",
            color="white", bbox=dict(boxstyle="round,pad=0.25", fc="#555", ec="none"))
    ax.text(0.03, 0.04, "cleaned seed", transform=ax.transAxes, fontsize=17, weight="bold",
            color="white", bbox=dict(boxstyle="round,pad=0.25", fc="#7a5230", ec="none"))
    ax.annotate("", xy=(0.26, 0.50), xytext=(0.12, 0.62), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#555", lw=2.5))
    ax.annotate("", xy=(0.24, 0.19), xytext=(0.17, 0.10), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#7a5230", lw=2.5))
    weigh = Image.open(out / "seed_weighing_lab.jpg").crop((450, 645, 1650, 1800))
    x0, y0, w, h = ax.get_position().bounds
    ins = fig.add_axes([x0 + w * 0.60, y0 + h * 0.47, w * 0.39, h * 0.52])
    ins.imshow(weigh); ins.set_xticks([]); ins.set_yticks([])
    for sp in ins.spines.values():
        sp.set_edgecolor("white"); sp.set_linewidth(4)
    ins.text(0.5, 0.03, "weighed (0.1 mg)", transform=ins.transAxes, ha="center", va="bottom",
             fontsize=15, weight="bold", color="white",
             bbox=dict(boxstyle="round,pad=0.25", fc="#333", ec="none"))
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
    legend = VPacker(align="center", pad=0, sep=5, children=[
        TextArea(f"Seed count = (weight − {SEED_W0 * 1e3:.2f} mg) ÷ {SEED_W1 * 1e3:.3f} mg",
                 textprops=dict(fontsize=15, weight="bold", color="#111")),
        TextArea(f"95 % interval ± {SEED_PI / SEED_W1:.0f} seeds",
                 textprops=dict(fontsize=13, color="#555"))])
    box = AnchoredOffsetbox(loc="upper right", child=legend, pad=0.5, borderpad=0,
                            frameon=True, bbox_to_anchor=(0.99, 0.44),
                            bbox_transform=ax.transAxes)
    box.patch.set(boxstyle="round,pad=0.3", fc="white", ec="#888", lw=1.5)
    ax.add_artist(box)


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
    res, freqs, genos = {}, {}, {}
    for pop, site in sites.items():
        sub = ind[(ind["EO"] == site) & (ind["n_functional_copies"] > 0)]
        g = np.stack(sub["genotype"].values)
        f = compute_local_f(np.stack(sub["fg_copies"].values))
        comp = p_compat_sporophytic_empirical(g, f, class_i, zyg, n_fathers=2000, rng=rng)
        res[pop] = dict(n=len(sub), compat=comp, alleles=int((f > 0).sum()))
        freqs[pop] = f
        genos[pop] = g                                   # 4 SRK copies per plant (tetraploid)
    order = np.argsort(-(freqs[42] + freqs[34]))
    order = [i for i in order if freqs[42][i] + freqs[34][i] > 0]
    letter = {i: chr(65 + k) for k, i in enumerate(order)}
    rows = [dict(allele=letter[i], srk_group=fg[i], **{f"freq_P{p}": float(freqs[p][i]) for p in sites})
            for i in order]
    tab = pd.DataFrame(rows)
    tab.to_csv(TABLES / "step32_genetic_example_P42_P34.tsv", sep="\t", index=False)

    # S-type colours as in Figure 1 (breeding_system_cartoon), A = commonest
    s_cols = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#17becf",
              "#bcbd22", "#aec7e8", "#ffbb78", "#98df8a", "#c5b0d5"]
    s_col = {i: s_cols[k % len(s_cols)] for k, i in enumerate(order)}
    col = {42: TYPE_COLORS["Mate-limited"], 34: TYPE_COLORS["Reproductively intact"]}
    name = {42: "P42 (EO70)\nmate-limited", 34: "P34 (EO18)\nintact"}
    # stage 1 = panel 1 headers only (S-types and plants per population),
    # 2 = + every plant's genotype and the share of the 3 commonest sets,
    # 3 = the whole figure (document); 1 – 2 are talk build-up frames with
    # the same layout, so nothing moves between clicks
    bbox = {}                                            # crop of the whole figure, for every stage

    def draw(stage: int, out: str) -> None:
        jitter = np.random.default_rng(2039)             # same dots in every frame
        fig = plt.figure(figsize=(19, 7.0))
        gs = fig.add_gridspec(1, 3, width_ratios=[2.6, 0.5, 0.85], wspace=0.3)

        # A — SRK genotype of every plant: 4 copies (tetraploid), colour = S-type;
        # plants grouped by their set of S-types, commonest set first
        ax = fig.add_subplot(gs[0])
        import collections
        ncol, x0s = 10, {42: 0.0, 34: 12.6}
        for pop in (42, 34):
            g = genos[pop]
            sets = [tuple(sorted({int(x) for x in r if x >= 0})) for r in g]
            cnt = collections.Counter(sets)
            idx = sorted(range(len(g)), key=lambda i: (-cnt[sets[i]], sets[i], tuple(sorted(g[i]))))
            for n, i in enumerate(idx if stage >= 2 else []):
                cx, cy = x0s[pop] + n % ncol + 0.5, -(n // ncol) - 0.5
                ax.add_patch(plt.Rectangle((cx - 0.46, cy - 0.46), 0.92, 0.92, fc="#F2F2F2",
                                           ec="#CFCFCF", lw=0.6))
                for k, a_ in enumerate(sorted(int(x) for x in g[i])):
                    if a_ < 0:
                        continue
                    ax.add_patch(plt.Circle((cx + (-0.21 if k % 2 == 0 else 0.21),
                                             cy + (0.21 if k < 2 else -0.21)), 0.18,
                                            fc=s_col.get(a_, "#999"), ec="white", lw=0.5))
            top3 = sum(c for _, c in cnt.most_common(3)) / len(g)
            f = freqs[pop]
            ax.text(x0s[pop] + 3.0, 3.85, name[pop].replace("\n", " — "), ha="left", va="top",
                    fontsize=15, color=col[pop])
            ax.text(x0s[pop] + 3.0, 2.75, f"{int((f > 0).sum())} S-types · {len(g)} plants",
                    ha="left", va="top", fontsize=15, weight="bold", color=col[pop])
            if stage >= 2:                               # under the plants: how repetitive
                rows = -(-len(g) // ncol)
                ax.text(x0s[pop], -rows - 0.35, f"3 S-type sets = {100 * top3:.0f} % of plants",
                        ha="left", va="top", fontsize=15, weight="bold", color=col[pop])
        ax.set_xlim(-0.1, 22.7); ax.set_ylim(-9.6, 5.4); ax.set_aspect("equal")   # room under the title
        ax.set_anchor("N")
        ax.axis("off")
        ax_geno = ax
        ax.set_title("1. S-locus genotyping (SRK): the S-types each plant carries", loc="left", weight="bold")

        # arrow: genotypes → mating (random mating within the population)
        # (a pollinator next to it: mating partners are drawn at random by pollinators)
        from restoration_cartoon import bee
        ax = fig.add_subplot(gs[1]); ax.axis("off"); ax.set_visible(stage >= 3)
        ax.set_xlim(0, 10); ax.set_ylim(0, 30); ax.set_aspect("equal"); ax.set_anchor("W")
        ax.annotate("", xy=(8.6, 15), xytext=(-3.5, 15), annotation_clip=False,
                    arrowprops=dict(arrowstyle="-|>,head_width=0.6,head_length=0.9", lw=4, color="#333"))
        bee(ax, 2.6, 19.2, s=1.3)
        ax.plot([-2.4, -0.6, 1.0], [17.0, 18.6, 18.9], ls=(0, (2, 2)), color="#777", lw=1.2)
        ax.text(2.6, 21.8, "random mating", ha="center", va="bottom", fontsize=15, weight="bold")
        ax.text(2.6, 13.2, "each plant meets\nfathers from its\npopulation; a shared\nS-type = rejected",
                ha="center", va="top", fontsize=11.5, color="#444", linespacing=1.15)

        # B — compatible mates per plant
        ax = fig.add_subplot(gs[2])
        x_mates = {42: 0.0, 34: 2.0}                              # spread apart to declutter
        for pop in (42, 34):
            xi = x_mates[pop]
            c = res[pop]["compat"]
            ax.scatter(xi + jitter.uniform(-0.32, 0.32, len(c)), c, s=34, color=col[pop], alpha=0.6,
                       linewidths=0)
            ax.plot([xi - 0.45, xi + 0.45], [c.mean()] * 2, color="black", lw=3)
            ax.text(xi + 0.5, c.mean(), f"{100 * c.mean():.0f} %", va="center", fontsize=15,
                    weight="bold")
        ax.set_xticks([x_mates[42], x_mates[34]])
        ax.set_xticklabels([name[p].replace(" (", "\n(") for p in (42, 34)], fontsize=16)
        for lab, pop in zip(ax.get_xticklabels(), (42, 34)):
            lab.set_color(col[pop])
        ax.tick_params(axis="y", labelsize=16, length=6, width=1.2)
        ax.tick_params(axis="x", length=6, width=1.2)
        ax.set_ylim(0, 1.40); ax.set_xlim(-0.55, 2.75)
        ax.set_yticks(np.arange(0, 1.01, 0.2))
        ax.spines["left"].set_bounds(0, 1)
        for pop in (42, 34):                                 # genetics agrees with the seed-set type
            ax.text(x_mates[pop], 1.115, "✓", ha="center", va="center", fontsize=26,
                    weight="bold", color=col[pop])
            ax.text(x_mates[pop], 1.045, "confirmed", ha="center", va="center", fontsize=12,
                    weight="bold", color=col[pop])
        ax_mates = ax
        from matplotlib.ticker import PercentFormatter
        ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
        ax.set_ylabel("Share of the population\na plant can mate with")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title("2. Compatible mates per plant", loc="left", weight="bold")
        ax.set_visible(stage >= 3)

        # boxed legend under the panels (as Figure 10)
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        h = [Patch(fc="#F2F2F2", ec="#CFCFCF", label="1 square = 1 genotyped plant;\n4 dots = its 4 SRK copies (tetraploid)"),
             Line2D([], [], marker="o", ls="", color=s_cols[0], ms=9,
                    label="color = S-type (as in Figure 1)"),
             Line2D([], [], marker="o", ls="", color="#888", ms=6, alpha=0.6,
                    label="dot = one plant; bar = population average")]
        h = h[:2] if stage == 2 else h if stage == 3 else []
        leg = fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=3,
                         fontsize=10.5, frameon=True, edgecolor="#999999", fancybox=False,
                         title="Legend", title_fontproperties={"weight": "bold", "size": 11.5},
                         columnspacing=1.8, handletextpad=0.5, borderpad=0.7)
        leg.get_frame().set_linewidth(1.0)
        leg.set_visible(bool(h))
        # type badges (as Figures 2 – 3) left of each population's header in panel A
        from type_icons import add_icon
        fig.canvas.draw()                                    # final layout of the equal-aspect panel
        to_fig = (ax_geno.transData + fig.transFigure.inverted()).transform
        w_over_h = fig.get_figwidth() / fig.get_figheight()
        for pop, kind in ((42, "Mate-limited"), (34, "Reproductively intact")):
            (xa, ya), (xb, yb) = to_fig([(x0s[pop] + 0.0, 0.25), (x0s[pop] + 2.7, 3.95)])
            side = min(xb - xa, (yb - ya) / w_over_h)
            add_icon(fig, kind, [xa, ya, side, side * w_over_h])
        # and above each population in panel B, with the genetic confirmation
        to_fig_b = (ax_mates.transData + fig.transFigure.inverted()).transform
        for pop, kind in ((42, "Mate-limited"), (34, "Reproductively intact")) if stage >= 3 else ():
            (cx, cy), = to_fig_b([(x_mates[pop], 1.28)])
            side = 0.046
            add_icon(fig, kind, [cx - side / 2, cy - side * w_over_h / 2, side, side * w_over_h])
        for ext in ("png", "pdf"):
            if stage == 3:
                bbox[ext] = fig.get_tightbbox(fig.canvas.get_renderer()).padded(0.1)
            fig.savefig(OUT / f"{out}.{ext}", dpi=220, bbox_inches=bbox[ext])
        plt.close(fig)

    draw(3, "genetic_example_P42_P34")                   # first: sets the crop
    draw(1, "genetic_example_P42_P34_step1")
    draw(2, "genetic_example_P42_P34_step2")
    print("[step32] genetic example:", {f"P{p}": dict(n=v["n"], alleles=v["alleles"],
          mean_compat=round(float(v["compat"].mean()), 3)) for p, v in res.items()})
    return tab


def materials_by_season() -> pd.DataFrame:
    """Materials shared by both parts of Seed_bulking (§ 5): per season,
    locations, EOs, slick spots (events) surveyed, reproductive plants
    counted, and plants measured for size and seed yield (LEPA DB; counts
    only, no coordinates). Demes and populations are results of Part 1,
    not materials, so they are not in this table."""
    import sqlite3
    from project_paths import LEPA_DB
    con = sqlite3.connect(f"file:{LEPA_DB}?mode=ro", uri=True)
    ev = pd.read_sql("""SELECT e.eventID, e.eventYear AS season, e.locationID, l.EOID,
                               e.organismQuantityFertile AS plants, e.populationID
                        FROM Events e JOIN Locations l USING(locationID)
                        WHERE e.populationID IS NOT NULL""", con)
    con.close()
    measured = load_occurrences()[["eventID", "year"]].merge(ev[["eventID"]], on="eventID")
    rows = []
    for season, g in list(ev.groupby("season")) + [("both", ev)]:
        n_meas = (measured["year"] == season).sum() if season != "both" else len(measured)
        rows.append(dict(season=season, EOs=g["EOID"].nunique(),
                         locations=g["locationID"].nunique(), slick_spots=g["eventID"].nunique(),
                         reproductive_plants=int(g["plants"].sum()),
                         plants_measured=int(n_meas)))
    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "step32_materials_by_season.tsv", sep="\t", index=False)
    return out


def main() -> None:
    OUT.mkdir(exist_ok=True)
    pops = population_table()
    pops.to_csv(TABLES / "step32_seed_bulking_populations.tsv", sep="\t", index=False)
    print(pops.groupby("type", observed=True)["populationID"].apply(list).to_string())
    m = plants_with_type(pops)
    plant_size_panel(m)
    seed_by_size(m, "seed_by_size_plants")
    seed_by_size(m, "seed_by_size_points", bands=False)       # animation step 1
    population_seed_yield(pops)
    population_seed_yield(pops, data=False, name="population_seed_yield_frame")   # animation step 1
    filter_funnel(pops)
    genetic_example(pops)
    print(seed_available(pops).to_string(index=False))
    print(materials_by_season().to_string(index=False))
    print(f"[step32] wrote tables/step32_seed_bulking_populations.tsv and figures in {OUT}/")


if __name__ == "__main__":
    main()
