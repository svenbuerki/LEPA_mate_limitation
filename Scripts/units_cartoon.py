"""units_cartoon.py — schematic of how the analysis units are built:
slick spot (event) → deme → population.

A cartoon with invented positions (no real site), drawn to the rules of
Steps 01 and 07:
  1. **Event** — one slick spot surveyed in one season; its flowering
     plants are counted.
  2. **Deme** — within a season, slick spots closer than 50 m (pollinator
     flight) are linked; each connected group is a deme.
  3. **Population** — slick spots of all seasons (so far 2025 + 2026)
     closer than 500 m are linked; each connected group is a population,
     inside which demes are redrawn each season.

Neutral colours, distinct from the BL and reproductive-type palettes, so
the figure fits the outline, the sampling design and the talk.
Output: figures/units_cartoon.{png,pdf}, with a copy in figures_seed_bulking/ for the talk
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


def hull(pts: np.ndarray, pad: float) -> np.ndarray:
    a = np.linspace(0, 2 * np.pi, 30, endpoint=False)
    cloud = (pts[:, None, :] + pad * np.stack([np.cos(a), np.sin(a)], -1)[None]).reshape(-1, 2)
    return cloud[ConvexHull(cloud).vertices]


def groups(pts: np.ndarray, d: float) -> np.ndarray:
    return connected_components(cdist(pts, pts) < d, directed=False)[1]


def draw_spots(ax, rng, plants=True):
    for (x, y), s, n in zip(CLOSE, SIZE, PLANTS):
        ax.add_patch(Ellipse((x, y), 2.4 * s, 1.6 * s, angle=rng.uniform(0, 180),
                             fc=SPOT, ec=SPOT_EDGE, lw=1.2, zorder=2))
        if plants:
            r = rng.uniform(0, 0.7, n) * s
            th = rng.uniform(0, 2 * np.pi, n)
            ax.scatter(x + r * np.cos(th), y + 0.6 * r * np.sin(th), s=14, color=PLANT, zorder=3)


def bee(ax, x, y, s=7.0, angle=20):
    """Small pollinator icon: striped body, two wings (data units = metres)."""
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
    rng = np.random.default_rng(7)
    fig, axes = plt.subplots(3, 1, figsize=(15, 19))        # one wide panel per row

    # 1 — events
    ax = axes[0]
    ax.set_xlim(0, 290); ax.set_ylim(0, 175)
    draw_spots(ax, np.random.default_rng(1))
    ax.annotate("one slick spot\n= one event", xy=CLOSE[0] + [0, 12], xytext=(18, 140),
                fontsize=13, arrowprops=dict(arrowstyle="->", color="#333"))
    ax.annotate(f"{PLANTS[4]} flowering plants\ncounted", xy=CLOSE[4] + [0, -10], xytext=(150, 2),
                fontsize=13, arrowprops=dict(arrowstyle="->", color="#333"))
    scale_bar(ax, 220, 160, 50, "50 m")
    ax.set_title("1. Slick spots (events) — 2025", loc="left", weight="bold", fontsize=18)
    ax.text(1.04, 0.5, "Each slick spot surveyed\nin a season is an event;\nits flowering plants\nare counted.",
            transform=ax.transAxes, va="center", fontsize=16)

    # 2 — demes
    ax = axes[1]
    ax.set_xlim(0, 290); ax.set_ylim(0, 175)
    lab = groups(CLOSE, 50)
    for k, g in enumerate(np.unique(lab)):
        pts = CLOSE[lab == g]
        h = hull(pts, 16)
        ax.fill(h[:, 0], h[:, 1], color=DEME_COLS[k % len(DEME_COLS)], alpha=0.55, zorder=1)
    d = cdist(CLOSE, CLOSE)
    for i in range(len(CLOSE)):
        for j in range(i + 1, len(CLOSE)):
            if d[i, j] < 50:
                ax.plot(*CLOSE[[i, j]].T, color=LINK, lw=1.5, ls="--", zorder=1.5)
    draw_spots(ax, np.random.default_rng(1), plants=False)
    from matplotlib.patches import FancyArrowPatch
    # pollinator flying between slick spots of one deme
    a, b = CLOSE[0], CLOSE[3]
    ax.add_patch(FancyArrowPatch(a + [6, 10], b + [-4, 10], connectionstyle="arc3,rad=-0.45",
                                 arrowstyle="-|>", mutation_scale=16, lw=1.6, color="#333",
                                 zorder=5))
    bee(ax, (a[0] + b[0]) / 2 + 2, 103, s=6.5, angle=-8)
    # too far: no flight between demes
    c, e = CLOSE[3], CLOSE[4]
    ax.add_patch(FancyArrowPatch(c + [8, -6], e + [-10, 2], connectionstyle="arc3,rad=0.25",
                                 arrowstyle="-|>", mutation_scale=14, lw=1.3, color="#C0392B",
                                 ls=(0, (4, 3)), zorder=5))
    mx, my = (c[0] + e[0]) / 2 + 4, (c[1] + e[1]) / 2 - 18
    ax.text(mx, my, "✕", color="#C0392B", fontsize=20, ha="center", va="center", zorder=6,
            weight="bold")
    ax.text(mx - 6, my + 30, "more than\n50 m apart", color="#C0392B", fontsize=11, ha="center",
            va="bottom")
    for g, name in zip(np.unique(lab), ("deme 1", "deme 2", "deme 3", "deme 4")):
        c = CLOSE[lab == g].mean(axis=0)
        ax.text(c[0], c[1] - 30 if name != "deme 3" else c[1] + 22, name, ha="center",
                fontsize=13, weight="bold")
    scale_bar(ax, 220, 160, 50, "50 m")
    ax.set_title("2. Demes — within a season (here 2025)", loc="left", weight="bold", fontsize=18)
    ax.text(1.04, 0.5, "Pollinators move between\nslick spots less than\n50 m apart, which then\nform one deme; a slick\n"
            "spot on its own is a\ndeme by itself.", transform=ax.transAxes, va="center",
            fontsize=16)

    # 3 — populations
    ax = axes[2]
    W, H = 4000, 2550
    ax.set_xlim(0, W); ax.set_ylim(0, H)
    r = np.random.default_rng(3)
    bx, by = 560, 780                                 # where the close-up of panels 1 – 2 sits
    close = CLOSE + [bx, by]                          # the same nine slick spots (2025), true scale
    pts, yr = [close], [2025] * len(CLOSE)
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
    is_close = np.arange(len(pts)) < len(CLOSE)
    plab = groups(pts, 500)
    for g in np.unique(plab):
        h = hull(pts[plab == g], 170)
        ax.fill(h[:, 0], h[:, 1], fc="#F4F1EA", ec="#333", lw=2, ls="--", zorder=1)
        c = pts[plab == g].mean(axis=0)
        ax.text(c[0], pts[plab == g][:, 1].max() + 230, f"population {g + 1}", ha="center",
                fontsize=13, weight="bold")
    m25, m26 = yr == 2025, yr == 2026
    far = ~is_close
    ax.scatter(*pts[far & m25].T, s=55, facecolors="white", edgecolors=SPOT_EDGE, linewidths=1.8,
               zorder=3, label="slick spot, 2025")
    ax.scatter(*pts[far & m26].T, s=55, color=SPOT_EDGE, zorder=3, label="slick spot, 2026")
    ax.scatter(*pts[is_close].T, s=14, facecolors="white", edgecolors=SPOT_EDGE, linewidths=1.2,
               zorder=3)                                  # the close-up's slick spots, smaller
    ax.add_patch(Rectangle((bx - 10, by - 10), 310, 195, fill=False, ec="#C0392B", lw=1.8,
                           zorder=4))
    ax.annotate("panels 1 – 2\n(the same nine\nslick spots)", xy=(bx + 150, by - 12),
                xytext=(bx + 30, by - 420), ha="center", va="top", color="#C0392B", fontsize=12,
                arrowprops=dict(arrowstyle="->", color="#C0392B"))
    ax.legend(loc="lower left", bbox_to_anchor=(1.02, 0.05), fontsize=14, frameon=False,
              handletextpad=0.3, markerscale=1.4)
    scale_bar(ax, W - 900, 120, 500, "500 m")
    ax.set_title("3. Populations — across seasons", loc="left", weight="bold", fontsize=18)
    ax.text(1.04, 0.62, "Slick spots of all seasons\n(here 2025 and 2026)\nless than 500 m apart\nform one population;\n"
            "its demes are redrawn\neach season.", transform=ax.transAxes, va="center", fontsize=16)

    for ax in axes:
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("#BBBBBB")
    fig.text(0.99, 0.005, "Schematic — positions are invented, not a real site.", ha="right",
             fontsize=12, color="#777", style="italic")
    fig.tight_layout(h_pad=3)
    talk = Path("figures_seed_bulking"); talk.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"units_cartoon.{ext}", dpi=220, bbox_inches="tight")
        fig.savefig(talk / f"units_cartoon.{ext}", dpi=220, bbox_inches="tight")   # talk copy
    plt.close(fig)
    print("[units_cartoon] wrote figures/units_cartoon.{png,pdf} (+ copy in figures_seed_bulking/)")


if __name__ == "__main__":
    main()
