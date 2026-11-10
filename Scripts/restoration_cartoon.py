"""restoration_cartoon.py — schematic of the restoration approach, with
slickspot peppergrass (Lepidium papilliferum) as the model.

Four steps, left to right:
  1. **Wild seed is scarce** — wild populations have limited fertility, so
     there is little seed to start a programme with.
  2. **Seed bulking** — the few wild seeds from source populations are
     grown and hand-crossed within each source, one row per mother, to
     multiply seed.
  3. **Informed breeding** — plants are genotyped at SRK; crosses join
     plants that share no S-allele, rare alleles first, so the seed carries
     many S-alleles in balanced numbers.
  4. **Restore & monitor** — that seed is sown into a mate-limited
     population: new S-alleles give each plant more compatible mates, seed
     set recovers; seed set and S-alleles are monitored.
A banner states what the whole programme rests on: rapidly identifying
candidate source populations — our approach, drawn in
source_screening_cartoon.py.

A cartoon in the style of units_cartoon.py (white panels, thin grey frame,
bold "N. Title" above, slick-spot soil colours, bee icon): invented plants
and alleles, icons rather than text. Reproductive-type colours as in the
seed-bulking talk; S-allele colours are their own palette.
Output: figures_seed_bulking/restoration_cartoon.{png,pdf}, plus
restoration_cartoon_step1–3.png (same size, later steps blank) to build a slide.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.transforms as mt
import numpy as np
from matplotlib.path import Path as MPath
from matplotlib.patches import (Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon,
                                Rectangle)

OUT = Path("figures_seed_bulking")
INTACT, MATE_LIM, OVER = "#2E7D32", "#F57C00", "#C62828"   # traffic light: green, orange, red
STEM, LEAF, FLOWER, FRUIT = "#3B6E3B", "#5E8C3A", "#FFFFFF", "#A0522D"   # fruit: sienna, as the seeds
FRUIT_EDGE = "#5D3A1A"
SPOT, SPOT_EDGE = "#E8D3A8", "#B08D57"
POT, POT_EDGE, SEED = "#C9A27E", "#7A5B3E", "#6D4C41"
NO = "#C62828"
ALLELE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
          "#8c564b", "#17becf", "#bcbd22"]           # S1 … S8
plt.rcParams.update({"font.size": 14})
RNG = np.random.default_rng(7)


# ---- icons (axes in data units 0–100) -------------------------------------
def lepa(ax, x, y, h=12, alleles=None, faded=False, fruits=True):
    """Slickspot peppergrass: short, much-branched stem, clusters of small
    white flowers, round fruits (True = on every branch, an int = on that
    many branches); optional S-allele dots at the base."""
    n_fruit = 5 if fruits is True else int(fruits)
    a = 0.35 if faded else 1.0
    ax.plot([x, x], [y, y + 0.45 * h], color=STEM, lw=2.2, alpha=a, zorder=3)
    for b, ang in enumerate((-50, -25, 0, 25, 50)):
        r = np.deg2rad(90 + ang)
        bx, by = x, y + 0.45 * h
        tx, ty = bx + 0.55 * h * np.cos(r), by + 0.55 * h * np.sin(r)
        ax.plot([bx, tx], [by, ty], color=STEM, lw=1.4, alpha=a, zorder=3)
        for _ in range(4):
            ax.add_patch(Circle((tx + RNG.normal(0, 0.06 * h), ty + RNG.normal(0, 0.05 * h)),
                                0.075 * h, fc=FLOWER, ec="#999", lw=0.5, alpha=a, zorder=4))
        if b in (2, 0, 4, 1, 3)[:n_fruit]:                # centre branch first
            ax.add_patch(Circle((tx - 0.08 * h, ty - 0.14 * h), 0.06 * h, fc=FRUIT, ec=FRUIT_EDGE,
                                lw=0.5, zorder=5))   # fruits stay visible on faded plants
    for s in (-1, 1):
        ax.add_patch(Polygon([[x, y], [x + s * 0.35 * h, y + 0.08 * h],
                              [x + s * 0.2 * h, y + 0.18 * h]], fc=LEAF, ec="none", alpha=a,
                             zorder=3))
    for k, al in enumerate(alleles or []):
        ax.add_patch(Circle((x + (k - 0.5) * 0.32 * h, y - 0.18 * h), 0.13 * h, fc=ALLELE[al],
                            ec="k", lw=0.6, alpha=a, zorder=5))


def seed_bag(ax, x, y, w, h, n, seed_cols=None, fc="#F5ECD9"):
    ax.add_patch(FancyBboxPatch((x - w / 2, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                fc=fc, ec="#7A6A55", lw=1.2, zorder=3))
    ax.plot([x - w / 2 - 0.2, x + w / 2 + 0.2], [y + h * 0.82] * 2, color="#7A6A55", lw=1.2,
            zorder=4)
    cols = seed_cols if seed_cols is not None else [SEED]
    sz = min(1.15, 0.16 * w)                             # seed length, data units
    pts = []                                             # spaced so seeds do not overlap
    for _ in range(400):
        if len(pts) == n:
            break
        q = (RNG.uniform(x - w / 2 + sz, x + w / 2 - sz), RNG.uniform(y + 0.8 * sz, y + h * 0.72))
        if all(np.hypot(q[0] - u, q[1] - v) > 1.05 * sz for u, v in pts):
            pts.append(q)
    for i, (sx, sy) in enumerate(pts):
        seed(ax, sx, sy, sz, RNG.uniform(0, 180), cols[i % len(cols)])


def seed(ax, x, y, size, angle, col=SEED):
    """One seed: an oval with a darker rim and a small notch (the hilum)."""
    ax.add_patch(Ellipse((x, y), size, 0.62 * size, angle=angle, fc=col, ec="#3E2723", lw=0.5,
                         zorder=4))
    a = np.deg2rad(angle)
    ax.add_patch(Circle((x + 0.42 * size * np.cos(a), y + 0.42 * size * np.sin(a)),
                        0.09 * size, fc="#3E2723", ec="none", zorder=5))


def arrow(ax, p, q, col="#444", lw=2, rad=0.0, ls="-", head="-|>"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=head, mutation_scale=16, color=col, lw=lw,
                                 connectionstyle=f"arc3,rad={rad}", ls=ls, zorder=6))


def bee(ax, x, y, s=1.6):
    for dx in (-0.35, 0.35):
        ax.add_patch(Ellipse((x + dx * s, y + 0.75 * s), 0.9 * s, 1.3 * s, fc="#DDEEFF",
                             ec="#555", lw=0.6, zorder=7))
    ax.add_patch(Ellipse((x, y), 2.2 * s, 1.2 * s, fc="#F4C430", ec="#222", lw=0.8, zorder=8))
    for dx in (-0.3, 0.25):
        ax.add_patch(Rectangle((x + dx * s - 0.12 * s, y - 0.55 * s), 0.24 * s, 1.1 * s,
                               fc="#222", zorder=9))
    ax.add_patch(Circle((x + 1.15 * s, y), 0.38 * s, fc="#222", zorder=9))


def brush(ax, x, y, s=1.0, angle=35):
    """Pollination brush: wooden handle, ferrule, yellow pollen tip."""
    tr = mt.Affine2D().rotate_deg_around(x, y, angle) + ax.transData
    ax.add_patch(Rectangle((x - 0.5 * s, y), 1.0 * s, 7 * s, fc="#B5835A", ec="#6D4C41", lw=0.6,
                           transform=tr, zorder=8))
    ax.add_patch(Rectangle((x - 0.6 * s, y - 1.2 * s), 1.2 * s, 1.4 * s, fc="#B0BEC5", ec="#555",
                           lw=0.6, transform=tr, zorder=8))
    ax.add_patch(Polygon([[x - 0.7 * s, y - 1.2 * s], [x + 0.7 * s, y - 1.2 * s],
                          [x, y - 3.4 * s]], fc="#F4C430", ec="#8D6E00", lw=0.6, transform=tr,
                         zorder=8))


def greenhouse(ax, x0, y0, w, h):
    """Glass house: walls, gable roof, glazing bars."""
    roof = 0.22 * h
    eave = y0 + h - roof
    ax.add_patch(Rectangle((x0, y0), w, h - roof, fc="#EEF5EA", ec="#7FA36B", lw=1.6, zorder=1))
    ax.add_patch(Polygon([[x0 - 1, eave], [x0 + w / 2, y0 + h], [x0 + w + 1, eave]],
                         fc="#E3EFF6", ec="#7FA36B", lw=1.6, zorder=1))
    for k in range(1, 4):
        xx = x0 + k * w / 4
        top = eave + roof * (1 - abs(xx - (x0 + w / 2)) / (w / 2))
        ax.plot([xx, xx], [eave, top], color="#A9C6A0", lw=1, zorder=1)


def tube(ax, x, y, s=1.0, fill="#90CAF9"):
    ax.add_patch(FancyBboxPatch((x - 1.3 * s, y), 2.6 * s, 8 * s,
                                boxstyle=f"round,pad=0,rounding_size={1.2 * s}", fc="white",
                                ec="#555", lw=1, zorder=5))
    ax.add_patch(FancyBboxPatch((x - 1.3 * s, y), 2.6 * s, 3.5 * s,
                                boxstyle=f"round,pad=0,rounding_size={1.2 * s}", fc=fill,
                                ec="none", zorder=5))
    ax.add_patch(Rectangle((x - 1.7 * s, y + 8 * s), 3.4 * s, 1.4 * s, fc="#E57373", ec="#555",
                           lw=0.8, zorder=6))


def dna(ax, x, y, h=12, w=3.2):
    t = np.linspace(0, 2.2 * np.pi, 80)
    ax.plot(x + w * np.sin(t), y + t / t[-1] * h, color="#1565C0", lw=2, zorder=5)
    ax.plot(x - w * np.sin(t), y + t / t[-1] * h, color="#EF6C00", lw=2, zorder=5)
    for tt in np.linspace(0.3, 2.0 * np.pi, 7):
        ax.plot([x + w * np.sin(tt), x - w * np.sin(tt)], [y + tt / t[-1] * h] * 2, color="#999",
                lw=1, zorder=4)


def clipboard(ax, x, y, s=1.0):
    """Monitoring: clipboard with a rising bar chart."""
    ax.add_patch(FancyBboxPatch((x - 4 * s, y), 8 * s, 10 * s, boxstyle="round,pad=0.2",
                                fc="#D7CCC8", ec="#6D4C41", lw=1, zorder=5))
    ax.add_patch(Rectangle((x - 3.2 * s, y + 0.8 * s), 6.4 * s, 8 * s, fc="white", ec="none",
                           zorder=5))
    ax.add_patch(Rectangle((x - 1.6 * s, y + 9.2 * s), 3.2 * s, 1.5 * s, fc="#90A4AE", ec="#555",
                           lw=0.6, zorder=6))
    for k in range(4):
        ax.add_patch(Rectangle((x - 2.4 * s + 1.6 * k * s, y + 1.6 * s), 1.1 * s,
                               (1.4 + 1.5 * k) * s, fc=INTACT, ec="none", zorder=6))


def cross_mark(ax, x, y, s=2.2, col=NO, lw=3):
    ax.plot([x - s, x + s], [y - s, y + s], color=col, lw=lw, zorder=9)
    ax.plot([x - s, x + s], [y + s, y - s], color=col, lw=lw, zorder=9)


def times_badge(ax, x, y, text):
    ax.add_patch(Circle((x, y), 4.2, fc="white", ec="#444", lw=1.4, zorder=7))
    ax.text(x, y, text, ha="center", va="center", fontsize=13, fontweight="bold", zorder=8)


def frame(ax, n, title, sub):
    ax.set_xlim(0, 100), ax.set_ylim(0, 100), ax.set_aspect("equal")
    ax.set_xticks([]), ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor("#BBBBBB")
    ax.set_title(f"{n}. {title}", loc="left", fontsize=18, fontweight="bold", pad=10)
    ax.text(50, 3.5, sub, ha="center", va="bottom", fontsize=13, color="#222")


# ---- panel 1: choose sources ------------------------------------------------
def panel_wild(ax):
    frame(ax, 1, "Wild seed is scarce", "limited fertility → little seed to start with")
    spots = [(24, 74, 30, 15, 3), (74, 77, 22, 12, 2), (30, 44, 24, 13, 2), (75, 46, 30, 15, 3)]
    for k, (x, y, w, h, n) in enumerate(spots):
        ax.add_patch(Ellipse((x, y), w, h, fc=SPOT, ec=SPOT_EDGE, lw=1.4, zorder=1))
        for dx in np.linspace(-(n - 1) * 4, (n - 1) * 4, n):
            lepa(ax, x + dx, y - 5, h=8, faded=dx != 0, fruits=2 if dx == 0 else 1)  # scarce, never none
        seed_bag(ax, x + w / 2 + 2, y - 7, 4.5, 5, 1 + k % 2)
    for (x1, y1), (x2, y2) in (((40, 74), (62, 77)), ((44, 44), (59, 46)), ((26, 64), (29, 53))):
        ax.plot([x1, x2], [y1, y2], color="#999", lw=1.2, ls=(0, (2, 3)), zorder=0)
    ax.text(51, 79.5, "isolated", ha="center", fontsize=10.5, color="#777", style="italic")
    arrow(ax, (50, 33), (50, 26), lw=2.2)
    seed_bag(ax, 50, 12, 9, 9, 4)
    ax.text(58, 16, "few\nseeds", ha="left", va="center", fontsize=12, fontweight="bold")


# ---- panel 2: seed bulking ---------------------------------------------
def panel_bulking(ax):
    frame(ax, 2, "Seed bulking", "grow, hand-cross within source, one row per mother")
    seed_bag(ax, 12, 62, 8, 8, 5)
    ax.text(12, 58, "few", ha="center", va="top", fontsize=12)
    arrow(ax, (18, 66), (26, 66))
    greenhouse(ax, 28, 36, 66, 52)
    for r, yy in enumerate((61, 49, 37)):
        ax.text(31, yy + 3.5, f"M{r + 1}", fontsize=11, color="#555", va="center")
        for k in range(5):
            xx = 42 + 11.5 * k
            ax.add_patch(Rectangle((xx - 3, yy), 6, 2.5, fc=POT, ec=POT_EDGE, lw=0.8, zorder=2))
            lepa(ax, xx, yy + 2.5, h=7.5)
    brush(ax, 70, 73, s=0.9)
    arrow(ax, (61, 34), (61, 27))
    for x in (34, 45, 56, 67, 78, 89):
        seed_bag(ax, x, 13, 7, 8, 12)
    times_badge(ax, 14, 22, "×10")


# ---- panel 3: informed breeding -------------------------------------------
PLANTS = [(1, 4), (2, 5), (1, 2), (0, 7), (3, 7), (5, 6)]


def panel_breeding(ax):
    frame(ax, 3, "Informed breeding", "genotype SRK, cross plants sharing no S-allele")
    pos = [(22, 58), (42, 58), (62, 58), (82, 58), (32, 33), (72, 33)]
    for (x, y), al in zip(pos, PLANTS):
        lepa(ax, x, y, h=11, alleles=list(al))
    ax.add_patch(Ellipse((10, 86), 5, 2.6, angle=30, fc=LEAF, ec=STEM, lw=0.8))  # leaf → tube → DNA
    arrow(ax, (14, 86), (19, 86), lw=1.5)
    tube(ax, 23, 81, s=0.8)
    arrow(ax, (27, 86), (32, 86), lw=1.5)
    dna(ax, 37, 80, h=12, w=2.6)
    for k in range(8):                                   # S-allele key
        ax.add_patch(Circle((48 + 6 * k, 89), 2.1, fc=ALLELE[k], ec="k", lw=0.6))
        ax.text(48 + 6 * k, 85.5, f"S{k + 1}", ha="center", va="top", fontsize=9)
    for i, j in [(0, 1), (2, 3), (4, 5), (1, 4), (3, 5)]:     # compatible crosses
        (x1, y1), (x2, y2) = pos[i], pos[j]
        ax.add_patch(FancyArrowPatch((x1, y1 + 10), (x2, y2 + 10), arrowstyle="-", color=INTACT,
                                     lw=2.4, connectionstyle="arc3,rad=0.25", zorder=2))
    (x1, y1), (x2, y2) = pos[0], pos[2]                  # share S2 → no cross
    ax.add_patch(FancyArrowPatch((x1, y1 + 12), (x2, y2 + 12), arrowstyle="-", color=NO, lw=2,
                                 ls=(0, (3, 2)), connectionstyle="arc3,rad=-0.2", zorder=2))
    cross_mark(ax, (x1 + x2) / 2, y1 + 16, s=1.8)
    brush(ax, 52, 40, s=0.9, angle=-30)
    arrow(ax, (85, 26), (85, 19))
    seed_bag(ax, 85, 9.5, 9, 8, 24, seed_cols=ALLELE)


# ---- panel 4: restore & monitor -------------------------------------------
def panel_restoration(ax):
    frame(ax, 4, "Restore & monitor", "new S-alleles → more mates → seed set")
    for x0, ok, al_sets in ((25, False, [(0, 1), (0, 1), (1, 0)]),
                            (73, True, [(0, 1), (2, 3), (4, 5), (6, 7), (0, 6)])):
        col = INTACT if ok else MATE_LIM
        ax.add_patch(Ellipse((x0, 50), 42, 34, fc=SPOT, ec=col, lw=3, zorder=1))
        ax.text(x0, 70, "after" if ok else "before", ha="center", fontsize=14, fontweight="bold",
                color=col)
        xs = np.linspace(x0 - 13, x0 + 13, len(al_sets))
        for k, (xx, al) in enumerate(zip(xs, al_sets)):
            lepa(ax, xx, 48 + (4 if k % 2 else -4), h=8, alleles=list(al), fruits=ok)
        if ok:                                           # bee flights between compatible plants
            for a, b in ((0, 1), (1, 2), (2, 3), (3, 4)):
                arrow(ax, (xs[a], 60), (xs[b], 60), col="#555", lw=1.2, rad=-0.35, ls="--",
                      head="-")
            bee(ax, x0 + 3, 63)
            seed_bag(ax, x0 + 14, 22, 7, 7, 14)
        else:
            bee(ax, x0, 63)
            seed_bag(ax, x0 + 12, 22, 7, 7, 2)
    arrow(ax, (46.5, 50), (51.5, 50), lw=2.5)
    seed_bag(ax, 49, 79, 7, 7, 14, seed_cols=ALLELE, fc="#E0F2F1")       # sown bred seed
    arrow(ax, (49, 78), (49, 69), col=INTACT)
    clipboard(ax, 10, 17, s=1.0)


PANELS = (panel_wild, panel_bulking, panel_breeding, panel_restoration)


def grid_figure(panels, n_shown: int, bottom: float = 0.02):
    """Four steps on two rows (1–2 above, 3–4 below), each panel drawn on a
    0–100 square; steps after n_shown left blank (same pixel size, so the
    frames build up one step per slide). Arrows 1 → 2, 2 → 3 (down to the
    next row) and 3 → 4. Shared with source_screening_cartoon.py."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 15.5 + 14 * bottom))
    flat = axes.flat
    for k, (ax, f) in enumerate(zip(flat, panels)):
        if k < n_shown:
            f(ax)
        else:
            ax.set_xlim(0, 100), ax.set_ylim(0, 100), ax.set_aspect("equal"), ax.axis("off")
            ax.set_title(" ", fontsize=18, pad=10)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=bottom, wspace=0.1, hspace=0.12)
    fig.canvas.draw()                                    # final layout before placing arrows
    to_fig = lambda ax, xy: fig.transFigure.inverted().transform(ax.transData.transform(xy))
    style = dict(transform=fig.transFigure, arrowstyle="-|>", mutation_scale=26, color="#333",
                 lw=3)
    for k, (a, b) in enumerate(((axes[0, 0], axes[0, 1]), (axes[1, 0], axes[1, 1]))):
        if 2 * k + 1 >= n_shown:                         # 1 → 2 and 3 → 4, side by side
            continue
        fig.add_artist(FancyArrowPatch(to_fig(a, (100, 50)), to_fig(b, (0, 50)), shrinkA=6,
                                       shrinkB=6, **style))
    if n_shown >= 3:                                     # 2 → 3: elbow connector
        (x0, y0), (x1, y1) = to_fig(axes[0, 1], (50, 0)), to_fig(axes[1, 0], (88, 100))
        y_mid = y0 - 0.012                               # just under the top row
        fig.add_artist(FancyArrowPatch(path=MPath([(x0, y0 - 0.003), (x0, y_mid), (x1, y_mid),
                                                  (x1, y1 + 0.004)]), **style))
    return fig


def draw(n_shown: int):
    global RNG
    RNG = np.random.default_rng(7)                       # same plants in every frame
    fig = grid_figure(PANELS, n_shown, bottom=0.07)
    if n_shown == len(PANELS):                          # what the programme rests on
        fig.text(0.5, 0.03, "★  The whole program rests on rapidly identifying candidate\n"
                 "source populations for seed bulking — our approach", ha="center", va="center",
                 fontsize=16, fontweight="bold", color="white", linespacing=1.3,
                 bbox=dict(boxstyle="round,pad=0.6", fc=INTACT, ec="none"))
    return fig


def main():
    OUT.mkdir(exist_ok=True)
    fig = draw(len(PANELS))
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"restoration_cartoon.{ext}", dpi=200 if ext == "png" else None)
    plt.close(fig)
    for n in range(1, len(PANELS)):                      # talk build: steps 1, 1–2, …
        fig = draw(n)
        fig.savefig(OUT / f"restoration_cartoon_step{n}.png", dpi=200)
        plt.close(fig)
    print(f"  wrote {OUT}/restoration_cartoon.png (+ _step1–3 frames for slides)")


if __name__ == "__main__":
    main()
