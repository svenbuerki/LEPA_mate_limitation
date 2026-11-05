"""units_cartoon.py — schematic of how the analysis units are built:
slick spot (event) → deme → population.

A cartoon with invented positions (no real site), drawn to the rules of
Steps 01 and 07, read top to bottom:
  1. **Slick spots (events)** — one frame split in two, 2025 | 2026, on the
     same ground: each slick spot surveyed in a season is an event and its
     flowering plants are counted; some slick spots are not above ground in
     a season, others appear.
  2. **Demes, per season** — within each season, slick spots closer than
     50 m (pollinator flight) are linked; each connected group is a deme.
     Demes are redrawn every season (here a new 2026 slick spot links two
     2025 demes).
  3. **Population, seasons combined** — slick spots of both seasons closer
     than 500 m are linked; each connected group is a population.

Neutral colours, distinct from the BL and reproductive-type palettes, so
the figure fits the outline, the sampling design and the talk.
Output: figures/units_cartoon.{png,pdf}, with a copy in figures_seed_bulking/ for the talk,
and a compact landscape version (the three steps left to right, 2025 above 2026) for slides
and Seed_bulking: figures_seed_bulking/units_cartoon_landscape.{png,pdf}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, Rectangle
from scipy.sparse.csgraph import connected_components
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist

OUT = Path("figures")
SPOT = "#E8D3A8"        # slick-spot soil
SPOT_EDGE = "#B08D57"
PLANT = "#3B6E3B"
DEME_COLS = ["#8FB8DE", "#F2B880", "#B5D99C", "#D3A4D9"]
LINK = "#555555"
plt.rcParams.update({"font.size": 14})

# ---- close-up (metres): 9 slick spots; within-50 m groups of 4, 3, 1 and 1
CLOSE = np.array([[30, 70], [62, 92], [70, 52], [100, 78],        # deme 1
                  [175, 40], [205, 62], [215, 28],                # deme 2
                  [150, 118],                                     # isolated
                  [255, 120]])                                    # isolated
SIZE = np.array([9, 6, 7, 5, 8, 5, 6, 4, 6])                      # half-width, m
PLANTS = np.array([14, 6, 9, 4, 11, 5, 7, 2, 6])
# 2026 on the same ground: slick spots 2 and 7 not above ground; two new
# ones — one between the 2025 demes 1 and 2 (links them), one near spot 9
KEEP26 = np.array([0, 2, 3, 4, 5, 7, 8])
NEW26 = np.array([[138, 60], [262, 150]])
CLOSE26 = np.vstack([CLOSE[KEEP26], NEW26])
SIZE26 = np.r_[SIZE[KEEP26], [5, 5]]
PLANTS26 = np.array([11, 7, 5, 9, 3, 3, 4, 5, 3])
GAP = 30                                                          # between the 2025 and 2026 halves
HALF = 290


def hull(pts: np.ndarray, pad: float) -> np.ndarray:
    a = np.linspace(0, 2 * np.pi, 30, endpoint=False)
    cloud = (pts[:, None, :] + pad * np.stack([np.cos(a), np.sin(a)], -1)[None]).reshape(-1, 2)
    return cloud[ConvexHull(cloud).vertices]


def groups(pts: np.ndarray, d: float) -> np.ndarray:
    return connected_components(cdist(pts, pts) < d, directed=False)[1]


def draw_spots(ax, rng, plants=True, pts=CLOSE, size=SIZE, n_pl=PLANTS, dx=0.0):
    for (x, y), sz, n in zip(pts, size, n_pl):
        x = x + dx
        ax.add_patch(Ellipse((x, y), 2.4 * sz, 1.6 * sz, angle=rng.uniform(0, 180),
                             fc=SPOT, ec=SPOT_EDGE, lw=1.2, zorder=2))
        if plants:
            r = rng.uniform(0, 0.7, n) * sz
            th = rng.uniform(0, 2 * np.pi, n)
            ax.scatter(x + r * np.cos(th), y + 0.6 * r * np.sin(th), s=14, color=PLANT, zorder=3)


def ghosts(ax, dx):
    """2025 slick spots not above ground in 2026: dashed outline."""
    for i in sorted(set(range(len(CLOSE))) - set(KEEP26)):
        (x, y), sz = CLOSE[i], SIZE[i]
        ax.add_patch(Ellipse((x + dx, y), 2.4 * sz, 1.6 * sz, fc="none", ec="#AAAAAA", lw=1.2,
                             ls=(0, (3, 2)), zorder=2))


def split_frame(ax):
    """One frame, two seasons: 2025 left, 2026 right (same ground)."""
    ax.set_xlim(0, 2 * HALF + GAP); ax.set_ylim(0, 185)
    ax.axvline(HALF + GAP / 2, color="#BBBBBB", lw=1.5)
    for x0, yr in ((HALF / 2, "2025"), (HALF + GAP + HALF / 2, "2026")):
        ax.text(x0, 176, yr, ha="center", va="top", fontsize=16, weight="bold", color="#555")


def draw_demes(ax, pts, dx, names):
    lab = groups(pts, 50)
    order = np.argsort([pts[lab == g][:, 0].min() for g in np.unique(lab)])
    for k, gi in enumerate(order):
        g = np.unique(lab)[gi]
        h = hull(pts[lab == g], 16)
        ax.fill(h[:, 0] + dx, h[:, 1], color=DEME_COLS[k % len(DEME_COLS)], alpha=0.55, zorder=1)
        c = pts[lab == g].mean(axis=0)
        top = pts[lab == g][:, 1].max()
        ax.text(c[0] + dx, top + 20 if top < 120 else pts[lab == g][:, 1].min() - 30, names[k],
                ha="center", fontsize=12.5, weight="bold")
    d = cdist(pts, pts)
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if d[i, j] < 50:
                ax.plot(pts[[i, j], 0] + dx, pts[[i, j], 1], color=LINK, lw=1.5, ls="--",
                        zorder=1.5)
    return len(np.unique(lab))


def bee(ax, x, y, s=7.0, angle=20):
    """Small pollinator icon: striped body, two wings (data units = meters)."""
    from matplotlib.patches import Ellipse as E
    import matplotlib.transforms as mt
    tr = mt.Affine2D().rotate_deg_around(x, y, angle) + ax.transData
    for dx in (-0.35, 0.35):
        ax.add_patch(E((x + dx * s, y + 0.75 * s), 0.9 * s, 1.3 * s, fc="#DDEEFF", ec="#555",
                       lw=0.8, alpha=0.95, transform=tr, zorder=6))
    ax.add_patch(E((x, y), 2.2 * s, 1.2 * s, fc="#F4C430", ec="#222", lw=1.0, transform=tr,
                   zorder=7))
    for dx in (-0.3, 0.25):
        ax.add_patch(Rectangle((x + dx * s - 0.12 * s, y - 0.55 * s), 0.24 * s, 1.1 * s,
                               fc="#222", transform=tr, zorder=8))
    ax.add_patch(plt.Circle((x + 1.15 * s, y), 0.38 * s, fc="#222", transform=tr, zorder=8))


def scale_bar(ax, x, y, length, label, fs=12):
    ax.plot([x, x + length], [y, y], color="black", lw=3)
    ax.text(x + length / 2, y + 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0]), label,
            ha="center", va="bottom", fontsize=fs)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    from matplotlib.patches import FancyArrowPatch
    fig, axes = plt.subplots(3, 1, figsize=(15, 21), gridspec_kw=dict(height_ratios=[1, 1, 1.55]))
    d26 = HALF + GAP                                      # x offset of the 2026 half

    # 1 — slick spots, both seasons
    ax = axes[0]
    split_frame(ax)
    draw_spots(ax, np.random.default_rng(1))
    draw_spots(ax, np.random.default_rng(2), pts=CLOSE26, size=SIZE26, n_pl=PLANTS26, dx=d26)
    ghosts(ax, d26)
    ax.annotate("one slick spot\n= one event", xy=CLOSE[0] + [0, 12], xytext=(8, 150),
                fontsize=12.5, arrowprops=dict(arrowstyle="->", color="#333"))
    ax.annotate(f"{PLANTS[4]} reproductive\nplants counted", xy=CLOSE[4] + [0, -10],
                xytext=(120, 4), fontsize=12.5, arrowprops=dict(arrowstyle="->", color="#333"))
    gx, gy = CLOSE[1] + [d26, 0]
    ax.annotate("not above\nground in 2026", xy=(gx, gy + 7), xytext=(gx - 40, 150),
                fontsize=11.5, color="#888", arrowprops=dict(arrowstyle="->", color="#888"))
    nx, ny = NEW26[0] + [d26, 0]
    ax.annotate("new in 2026", xy=(nx, ny - 7), xytext=(nx - 10, 6), fontsize=11.5,
                color="#555", ha="center", arrowprops=dict(arrowstyle="->", color="#555"))
    scale_bar(ax, 2 * HALF + GAP - 70, 8, 50, "50 m")
    ax.set_title("1. Survey slick spots holding reproductive plants — each season", loc="left",
                 weight="bold", fontsize=18)
    ax.text(1.03, 0.5, "Each slick spot surveyed\nin a season is an event;\nits reproductive plants\n"
            "are counted. Some slick\nspots are not above\nground every season.",
            transform=ax.transAxes, va="center", fontsize=15)

    # 2 — demes, per season
    ax = axes[1]
    split_frame(ax)
    draw_demes(ax, CLOSE, 0, ("deme 1", "deme 2", "deme 3", "deme 4"))
    draw_demes(ax, CLOSE26, d26, ("deme 1", "deme 2", "deme 3", "deme 4"))
    draw_spots(ax, np.random.default_rng(1), plants=False)
    draw_spots(ax, np.random.default_rng(2), plants=False, pts=CLOSE26, size=SIZE26,
               n_pl=PLANTS26, dx=d26)
    a, b = CLOSE[0], CLOSE[3]                             # pollinator within a 2025 deme
    ax.add_patch(FancyArrowPatch(a + [6, 10], b + [-4, 10], connectionstyle="arc3,rad=-0.45",
                                 arrowstyle="-|>", mutation_scale=16, lw=1.6, color="#333",
                                 zorder=5))
    bee(ax, (a[0] + b[0]) / 2 + 2, 103, s=6.5, angle=-8)
    c, e = CLOSE[3], CLOSE[4]                             # too far in 2025
    ax.add_patch(FancyArrowPatch(c + [8, -6], e + [-10, 2], connectionstyle="arc3,rad=0.25",
                                 arrowstyle="-|>", mutation_scale=14, lw=1.3, color="#C0392B",
                                 ls=(0, (4, 3)), zorder=5))
    ax.text((c[0] + e[0]) / 2 + 4, (c[1] + e[1]) / 2 - 18, "✕", color="#C0392B", fontsize=20,
            ha="center", va="center", zorder=6, weight="bold")
    ax.text((c[0] + e[0]) / 2 - 2, 4, "more than 50 m apart", color="#C0392B", fontsize=11,
            ha="center", va="bottom")
    ax.text(NEW26[0][0] + d26, 4, "a new slick spot links\ntwo 2025 demes", color="#555",
            fontsize=11, ha="center", va="bottom")
    ax.set_title("2. Determine the mating pool (demes) — each season", loc="left", weight="bold",
                 fontsize=18)
    ax.text(1.03, 0.5, "Plants less than 50 m\napart (pollinator flight)\ncan mate: their slick\n"
            "spots form one deme,\nthe mating pool, redrawn\neach season.", transform=ax.transAxes,
            va="center", fontsize=15)

    # 3 — populations
    ax = axes[2]
    W, H = 4000, 2550
    ax.set_xlim(0, W); ax.set_ylim(0, H)
    r = np.random.default_rng(3)
    bx, by = 560, 780                                 # where the close-up of panels 1 – 2 sits
    close = np.vstack([CLOSE, CLOSE26]) + [bx, by]    # panels 1 – 2, both seasons, true scale
    pts, yr = [close], [2025] * len(CLOSE) + [2026] * len(CLOSE26)
    # the rest of population 1 around the close-up (outside its box)
    extra = []
    while len(extra) < 7:
        p = r.normal([bx + 145, by + 90], 230)
        if not (bx - 40 < p[0] < bx + 330 and by - 40 < p[1] < by + 215):
            extra.append(p)
    pts.append(np.array(extra)); yr += list(r.choice([2025, 2026], size=len(extra)))
    clusters = [((1350, 1550), 7, 150), ((2900, 700), 10, 170), ((3100, 1900), 6, 140)]
    for (cx, cy), n, sd in clusters:
        p = r.normal([cx, cy], sd, size=(n, 2))
        pts.append(p); yr += list(r.choice([2025, 2026], size=n))
    pts = np.vstack(pts); yr = np.array(yr)
    is_close = np.arange(len(pts)) < len(CLOSE) + len(CLOSE26)
    plab = groups(pts, 500)
    for g in np.unique(plab):
        h = hull(pts[plab == g], 170)
        ax.fill(h[:, 0], h[:, 1], fc="#F4F1EA", ec="#333", lw=2, ls="--", zorder=1)
        c = pts[plab == g].mean(axis=0)
        ax.text(c[0], pts[plab == g][:, 1].max() + 230, f"population {g + 1}", ha="center",
                fontsize=13, weight="bold")
    # too far: population 1 (the close-up's) and population 4, across the open middle
    from matplotlib.patches import FancyArrowPatch
    g0 = plab[0]
    g1 = plab[np.argmin(cdist(pts, [[2900, 700]]))]
    dd = cdist(pts[plab == g0], pts[plab == g1])
    i0, j0 = np.unravel_index(np.argmin(dd), dd.shape)
    pa, pb = pts[plab == g0][i0], pts[plab == g1][j0]
    ax.add_patch(FancyArrowPatch(pa, pb, connectionstyle="arc3,rad=-0.2", arrowstyle="-|>",
                                 mutation_scale=18, lw=2, color="#C0392B", ls=(0, (4, 3)),
                                 shrinkA=12, shrinkB=12, zorder=5))
    mid = (pa + pb) / 2 + [0, 90]
    ax.text(*mid, "✕", color="#C0392B", fontsize=24, ha="center", va="center", zorder=6,
            weight="bold")
    ax.text(mid[0], mid[1] + 120, "more than 500 m apart", color="#C0392B", fontsize=13,
            ha="center", va="bottom", zorder=6)
    m25, m26 = yr == 2025, yr == 2026
    far = ~is_close
    ax.scatter(*pts[far & m25].T, s=55, facecolors="white", edgecolors=SPOT_EDGE, linewidths=1.8,
               zorder=3, label="slick spot, 2025")
    ax.scatter(*pts[far & m26].T, s=55, color=SPOT_EDGE, zorder=3, label="slick spot, 2026")
    ax.scatter(*pts[is_close & m25].T, s=14, facecolors="white", edgecolors=SPOT_EDGE,
               linewidths=1.2, zorder=3)                  # the close-up's slick spots, smaller
    ax.scatter(*pts[is_close & m26].T, s=10, color=SPOT_EDGE, zorder=3)
    ax.add_patch(Rectangle((bx - 10, by - 10), 310, 195, fill=False, ec="#C0392B", lw=1.8,
                           zorder=4))
    ax.annotate("panels 1 – 2\n(both seasons)", xy=(bx + 150, by - 12),
                xytext=(bx + 30, by - 420), ha="center", va="top", color="#C0392B", fontsize=12,
                arrowprops=dict(arrowstyle="->", color="#C0392B"))
    ax.legend(loc="lower left", bbox_to_anchor=(1.02, 0.05), fontsize=14, frameon=False,
              handletextpad=0.3, markerscale=1.4)
    scale_bar(ax, W - 900, 120, 500, "500 m")
    ax.set_title("3. Determine populations — the units to prioritize", loc="left", weight="bold",
                 fontsize=18)
    ax.text(1.03, 0.62, "Slick spots of all seasons\n(here 2025 and 2026)\nless than 500 m apart\n(maximum pollinator\n"
            "dispersal) form one\npopulation, the unit\nto prioritize.", transform=ax.transAxes, va="center", fontsize=15)

    for ax in axes:
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("#BBBBBB")
    fig.text(0.99, 0.005, "Schematic — positions are invented, not a real site.", ha="right",
             fontsize=12, color="#777", style="italic")
    fig.tight_layout(h_pad=3.2)
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    steps = [((axes[0], axes[1]), [HALF / 2, HALF + GAP + HALF / 2], "group slick spots < 50 m apart, each season"),
             ((axes[1], axes[2]), None, "combine both seasons; link slick spots\n< 500 m apart (maximum pollinator dispersal)")]
    for (a0, a1), xs, lab in steps:
        y0 = inv.transform(a0.transAxes.transform((0, 0)))[1] - 0.004
        y1 = inv.transform(a1.transAxes.transform((0, 1)))[1] + 0.022
        if xs is None:                                    # both halves merge (elbow connectors)
            from matplotlib.lines import Line2D
            xs_f = [inv.transform(a0.transData.transform((x, 0)))[0] for x in (HALF / 2, HALF + GAP + HALF / 2)]
            xm = inv.transform(a0.transData.transform((HALF + GAP / 2, 0)))[0]
            y_mid = (y0 + y1) / 2
            for xf in xs_f:                               # down from each season, across to the middle
                fig.add_artist(Line2D([xf, xf, xm], [y0, y_mid, y_mid], transform=fig.transFigure,
                                      lw=2.5, color="#333", solid_joinstyle="miter"))
            fig.add_artist(FancyArrowPatch((xm, y_mid), (xm, y1), transform=fig.transFigure,
                                           arrowstyle="-|>", mutation_scale=22, lw=2.5,
                                           color="#333", shrinkA=0))
            fig.text(xs_f[1] + 0.025, y_mid, lab, fontsize=12.5, style="italic", color="#333",
                     va="center")
        else:
            for x in xs:
                xf = inv.transform(a0.transData.transform((x, 0)))[0]
                fig.add_artist(FancyArrowPatch((xf, y0), (xf, y1), transform=fig.transFigure,
                                               arrowstyle="-|>", mutation_scale=22, lw=2.5,
                                               color="#333"))
            xm = inv.transform(a0.transAxes.transform((0.5, 0)))[0]
            fig.text(xm, (y0 + y1) / 2, lab, fontsize=12.5, style="italic",
                     color="#333", ha="center", va="center",
                     bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none"))
    talk = Path("figures_seed_bulking"); talk.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"units_cartoon.{ext}", dpi=220, bbox_inches="tight")
        fig.savefig(talk / f"units_cartoon.{ext}", dpi=220, bbox_inches="tight")   # talk copy
    plt.close(fig)
    print("[units_cartoon] wrote figures/units_cartoon.{png,pdf} (+ copy in figures_seed_bulking/)")


# ---- landscape version (talk): the three steps left to right ---------------------------
def year_frame(ax, yr):
    ax.set_xlim(0, HALF); ax.set_ylim(0, 185); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("#BBBBBB")
    ax.text(8, 176, yr, ha="left", va="top", fontsize=15, weight="bold", color="#555")


def landscape_spots(a25, a26):
    year_frame(a25, "2025"); year_frame(a26, "2026")
    draw_spots(a25, np.random.default_rng(1))
    draw_spots(a26, np.random.default_rng(2), pts=CLOSE26, size=SIZE26, n_pl=PLANTS26)
    ghosts(a26, 0)
    a25.annotate("one slick spot\n= one event", xy=CLOSE[0] + [0, -12], xytext=(4, 8),
                 fontsize=11.5, arrowprops=dict(arrowstyle="->", color="#333"))
    a25.annotate(f"{PLANTS[4]} reproductive\nplants counted", xy=CLOSE[4] + [0, 10],
                 xytext=(150, 150), fontsize=11.5, arrowprops=dict(arrowstyle="->", color="#333"))
    gx, gy = CLOSE[1]
    a26.annotate("not above\nground", xy=(gx, gy + 7), xytext=(gx - 25, 150), fontsize=11,
                 color="#888", arrowprops=dict(arrowstyle="->", color="#888"))
    nx, ny = NEW26[0]
    a26.annotate("new in 2026", xy=(nx, ny - 7), xytext=(nx - 10, 6), fontsize=11,
                 color="#555", ha="center", arrowprops=dict(arrowstyle="->", color="#555"))
    scale_bar(a26, HALF - 70, 8, 50, "50 m", fs=11)


def landscape_demes(a25, a26):
    from matplotlib.patches import FancyArrowPatch
    year_frame(a25, "2025"); year_frame(a26, "2026")
    draw_demes(a25, CLOSE, 0, ("deme 1", "deme 2", "deme 3", "deme 4"))
    draw_demes(a26, CLOSE26, 0, ("deme 1", "deme 2", "deme 3", "deme 4"))
    draw_spots(a25, np.random.default_rng(1), plants=False)
    draw_spots(a26, np.random.default_rng(2), plants=False, pts=CLOSE26, size=SIZE26,
               n_pl=PLANTS26)
    a, b = CLOSE[0], CLOSE[3]                             # pollinator within a 2025 deme
    a25.add_patch(FancyArrowPatch(a + [6, 10], b + [-4, 10], connectionstyle="arc3,rad=-0.45",
                                  arrowstyle="-|>", mutation_scale=16, lw=1.6, color="#333",
                                  zorder=5))
    bee(a25, (a[0] + b[0]) / 2 + 2, 103, s=6.5, angle=-8)
    c, e = CLOSE[3], CLOSE[4]                             # too far in 2025
    a25.add_patch(FancyArrowPatch(c + [8, -6], e + [-10, 2], connectionstyle="arc3,rad=0.25",
                                  arrowstyle="-|>", mutation_scale=14, lw=1.3, color="#C0392B",
                                  ls=(0, (4, 3)), zorder=5))
    a25.text((c[0] + e[0]) / 2 + 4, (c[1] + e[1]) / 2 - 18, "✕", color="#C0392B", fontsize=18,
             ha="center", va="center", zorder=6, weight="bold")
    a25.text((c[0] + e[0]) / 2 - 2, 4, "more than 50 m apart", color="#C0392B", fontsize=10.5,
             ha="center", va="bottom")
    a26.text(NEW26[0][0], 4, "a new slick spot links two 2025 demes", color="#555",
             fontsize=10.5, ha="center", va="bottom")


def landscape_populations(ax):
    from matplotlib.patches import FancyArrowPatch
    W3, H3 = 3000, 2900
    ax.set_xlim(0, W3); ax.set_ylim(0, H3); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("#BBBBBB")
    r = np.random.default_rng(4)
    bx, by = 520, 1150                                    # the close-up of panels 1 – 2
    close = np.vstack([CLOSE, CLOSE26]) + [bx, by]
    pts, yr = [close], [2025] * len(CLOSE) + [2026] * len(CLOSE26)
    extra = []
    while len(extra) < 7:
        q = r.normal([bx + 145, by + 90], 210)
        if not (bx - 40 < q[0] < bx + 330 and by - 40 < q[1] < by + 215) and q[0] > 60:
            extra.append(q)
    pts.append(np.array(extra)); yr += list(r.choice([2025, 2026], size=len(extra)))
    clusters = [((1050, 2250), 7, 140), ((2400, 1100), 10, 150), ((2450, 2200), 6, 120),
                ((1300, 300), 1, 0)]                      # the last: a one-slick-spot population
    for (cx, cy), n, sd in clusters:
        q = r.normal([cx, cy], sd, size=(n, 2)) if sd else np.array([[cx, cy]])
        pts.append(q); yr += list(r.choice([2025, 2026], size=n))
    pts = np.vstack(pts); yr = np.array(yr)
    is_close = np.arange(len(pts)) < len(CLOSE) + len(CLOSE26)
    plab = groups(pts, 500)
    for g in np.unique(plab):
        h = hull(pts[plab == g], 150)
        ax.fill(h[:, 0], h[:, 1], fc="#F4F1EA", ec="#333", lw=2, ls="--", zorder=1)
        c = pts[plab == g].mean(axis=0)
        ax.text(c[0], pts[plab == g][:, 1].max() + 200, f"population {g + 1}", ha="center",
                fontsize=12.5, weight="bold")
    g0 = plab[0]
    g1 = plab[np.argmin(cdist(pts, [[2400, 1100]]))]
    dd = cdist(pts[plab == g0], pts[plab == g1])
    i0, j0 = np.unravel_index(np.argmin(dd), dd.shape)
    pa, pb = pts[plab == g0][i0], pts[plab == g1][j0]
    ax.add_patch(FancyArrowPatch(pa, pb, connectionstyle="arc3,rad=-0.2", arrowstyle="-|>",
                                 mutation_scale=18, lw=2, color="#C0392B", ls=(0, (4, 3)),
                                 shrinkA=12, shrinkB=12, zorder=5))
    mid = (pa + pb) / 2 + [0, 80]
    ax.text(*mid, "✕", color="#C0392B", fontsize=22, ha="center", va="center", zorder=6,
            weight="bold")
    ax.text(mid[0], mid[1] + 110, "more than 500 m apart", color="#C0392B", fontsize=12,
            ha="center", va="bottom", zorder=6)
    m25, m26, far = yr == 2025, yr == 2026, ~is_close
    ax.scatter(*pts[far & m25].T, s=50, facecolors="white", edgecolors=SPOT_EDGE, linewidths=1.8,
               zorder=3, label="slick spot, 2025")
    ax.scatter(*pts[far & m26].T, s=50, color=SPOT_EDGE, zorder=3, label="slick spot, 2026")
    ax.scatter(*pts[is_close & m25].T, s=12, facecolors="white", edgecolors=SPOT_EDGE,
               linewidths=1.1, zorder=3)
    ax.scatter(*pts[is_close & m26].T, s=9, color=SPOT_EDGE, zorder=3)
    ax.add_patch(Rectangle((bx - 10, by - 10), 310, 195, fill=False, ec="#C0392B", lw=1.8,
                           zorder=4))
    ax.annotate("panels 1 – 2", xy=(bx + 150, by - 12), xytext=(bx + 250, by - 450),
                ha="center", va="top", color="#C0392B", fontsize=11.5,
                arrowprops=dict(arrowstyle="->", color="#C0392B"))
    ax.legend(loc="lower right", fontsize=12, frameon=False, handletextpad=0.3,
              markerscale=1.2, borderaxespad=0.3, bbox_to_anchor=(1.0, 0.07))
    scale_bar(ax, W3 - 700, 90, 500, "500 m", fs=11)


def main_landscape() -> None:
    """Landscape, compact: slick spots → demes (2025 above 2026) → populations."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyArrowPatch
    cw, fh = 5.0, 5.0 * 185 / HALF                        # column width, year-frame height
    m, title, cap = 0.3, 0.95, 1.55
    h_col = 2 * fh + 0.25
    g1, g2, w3 = 1.25, 1.45, h_col * 3000 / 2900           # panel 3 as tall as the columns
    W, H = 2 * m + 2 * cw + g1 + g2 + w3, m + title + h_col + cap + 0.3
    fig = plt.figure(figsize=(W, H))
    def box(x, y_top, w, h):
        return fig.add_axes([x / W, 1 - (y_top + h) / H, w / W, h / H])
    y0 = m + title
    x1, x2, x3 = m, m + cw + g1, m + 2 * cw + g1 + g2
    p1 = (box(x1, y0, cw, fh), box(x1, y0 + fh + 0.25, cw, fh))
    p2 = (box(x2, y0, cw, fh), box(x2, y0 + fh + 0.25, cw, fh))
    p3 = box(x3, y0 + (h_col - w3 * 2900 / 3000) / 2, w3, w3 * 2900 / 3000)
    landscape_spots(*p1); landscape_demes(*p2); landscape_populations(p3)
    heads = ((x1, "1. Survey slick spots holding\n    reproductive plants"),
             (x2, "2. Determine the mating pool\n    (demes)"),
             (x3, "3. Determine populations —\n    the units to prioritize"))
    for x, s in heads:
        fig.text(x / W, 1 - (m + 0.8) / H, s, ha="left", va="baseline", fontsize=18,
                 weight="bold", linespacing=1.15)
    caps = ((x1, cw, "field survey, two seasons",
             "Each slick spot surveyed in a season is an event;\nits reproductive plants "
             "are counted. Some slick spots\nare not above ground every season."),
            (x2, cw, "pollinator flight, < 50 m",
             "Plants less than 50 m apart can mate: their\nslick spots form one deme, the "
             "mating pool,\nredrawn each season."),
            (x3, w3, "maximum pollinator dispersal, < 500 m",
             "Slick spots of both seasons within the maximum\npollinator dispersal (500 m) form "
             "one population —\nthe unit for seed collection, restoration and prioritization."))
    for x, w, tech, s in caps:
        xc = (x + w / 2) / W
        fig.text(xc, 1 - (y0 + h_col + 0.32) / H, f"technique:  {tech}", ha="center",
                 va="center", fontsize=13, color="white", weight="bold",
                 bbox=dict(boxstyle="round,pad=0.4", fc="#546E7A", ec="none"))
        fig.text(xc, 1 - (y0 + h_col + 0.7) / H, s, ha="center", va="top", fontsize=13,
                 linespacing=1.25)
    fig.canvas.draw()
    to_f = lambda ax, xy: fig.transFigure.inverted().transform(ax.transAxes.transform(xy))
    st = dict(transform=fig.transFigure, arrowstyle="-|>", mutation_scale=22, lw=2.5,
              color="#333")
    for a, b in zip(p1, p2):                              # 1 → 2, one arrow per season
        fig.add_artist(FancyArrowPatch(to_f(a, (1, 0.5)), to_f(b, (0, 0.5)), shrinkA=4,
                                       shrinkB=4, **st))
    xa, ya = to_f(p1[0], (1, 0.5)); xb, _ = to_f(p2[0], (0, 0.5))
    fig.text((xa + xb) / 2, ya + 0.03, "< 50 m", ha="center", va="bottom", fontsize=12.5,
             style="italic")
    (xr, y25), (_, y26) = to_f(p2[0], (1, 0.5)), to_f(p2[1], (1, 0.5))   # 2 → 3, merged
    xm = xr + 0.35 * (g2 / W)
    for y in (y25, y26):
        fig.add_artist(Line2D([xr + 0.003, xm, xm], [y, y, (y25 + y26) / 2],
                              transform=fig.transFigure, lw=2.5, color="#333",
                              solid_joinstyle="miter"))
    x3l, _ = to_f(p3, (0, 0.5))
    fig.add_artist(FancyArrowPatch((xm, (y25 + y26) / 2), (x3l, (y25 + y26) / 2), shrinkA=0,
                                   shrinkB=4, **st))
    fig.text((xm + x3l) / 2 + 0.002, (y25 + y26) / 2 - 0.03, "< 500 m,\nboth\nseasons",
             ha="center", va="top", fontsize=12.5, style="italic", linespacing=1.1)
    fig.text(1 - m / W, 0.012, "Schematic — positions are invented, not a real site.",
             ha="right", fontsize=11.5, color="#777", style="italic")
    talk = Path("figures_seed_bulking"); talk.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(talk / f"units_cartoon_landscape.{ext}", dpi=220)
    plt.close(fig)
    print("[units_cartoon] wrote figures_seed_bulking/units_cartoon_landscape.{png,pdf}")


if __name__ == "__main__":
    main()
    main_landscape()
