"""source_screening_cartoon.py — schematic of our approach to rapidly
identify candidate source populations for seed bulking, with slickspot
peppergrass (Lepidium papilliferum) as the model. Companion of
restoration_cartoon.py, whose programme rests on this step.

Four steps, left to right:
  1. **Measure every plant** — crown size and seeds of every flowering
     plant, in two seasons: fast and low cost, no genotyping.
  2. **Seeds expected for size** — all plants together give the seed a
     plant of a given size is expected to set (band: half to twice).
  3. **Compare each population** — each population's seed set relative to
     its plants' size, per season: below (mate-limited), as expected
     (reproductively intact) or above (over-producing).
  4. **Candidate sources** — intact populations go to seed bulking;
     mate-limited ones are restoration targets; over-producing ones are
     excluded; SRK genotypes are the next check.

Same icons and colours as restoration_cartoon.py; invented plants and
values. Output: figures_seed_bulking/source_screening_cartoon.{png,pdf},
plus source_screening_cartoon_step1–3.png (same size, later steps blank).
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, Rectangle

import restoration_cartoon as rc
from restoration_cartoon import (INTACT, MATE_LIM, OVER, OUT, SPOT, SPOT_EDGE, arrow,
                                 cross_mark, dna, frame, greenhouse, lepa, seed_bag)

GREY_PT = "#9E9E9E"


def stopwatch(ax, x, y, r=4.5):
    ax.add_patch(Rectangle((x - 0.9, y + r), 1.8, 1.4, fc="#555", zorder=6))
    ax.add_patch(Circle((x, y), r, fc="white", ec="#333", lw=1.8, zorder=6))
    for a in np.deg2rad(np.arange(0, 360, 30)):
        ax.plot([x + 0.8 * r * np.cos(a), x + 0.95 * r * np.cos(a)],
                [y + 0.8 * r * np.sin(a), y + 0.95 * r * np.sin(a)], color="#333", lw=0.8, zorder=7)
    ax.plot([x, x + 0.55 * r * np.cos(np.deg2rad(40))], [y, y + 0.55 * r * np.sin(np.deg2rad(40))],
            color=INTACT, lw=2.2, zorder=7)
    ax.plot([x, x], [y, y + 0.7 * r], color="#333", lw=1.6, zorder=7)


def panel_measure(ax):
    frame(ax, 1, "Measure every plant", "crown size + seeds, every plant, two seasons")
    ax.add_patch(Ellipse((50, 74), 70, 22, fc=SPOT, ec=SPOT_EDGE, lw=1.4, zorder=1))
    for k, (dx, h) in enumerate(((-26, 7), (-14, 10), (-2, 6), (10, 12), (22, 8))):
        lepa(ax, 50 + dx, 67, h=h)
    ax.text(50, 89, "every flowering plant", ha="center", fontsize=11.5, color="#555",
            style="italic")
    # one plant on the measuring board: crown width and height
    ax.add_patch(Rectangle((8, 25), 34, 2, fc="#ECEFF1", ec="#90A4AE", lw=0.8, zorder=2))
    ax.add_patch(Rectangle((8, 25), 2, 30, fc="#ECEFF1", ec="#90A4AE", lw=0.8, zorder=2))
    for t in range(10):
        ax.plot([10, 11.2], [28 + 2.8 * t] * 2, color="#607D8B", lw=0.8, zorder=3)
        ax.plot([11.5 + 3 * t] * 2, [27, 28.2], color="#607D8B", lw=0.8, zorder=3)
    lepa(ax, 26, 27, h=20)
    ax.add_patch(FancyArrowPatch((15, 51), (37, 51), arrowstyle="<|-|>", mutation_scale=10,
                                 color="#333", lw=1.2, zorder=6))
    ax.text(26, 52.5, "width", ha="center", va="bottom", fontsize=10.5)
    ax.add_patch(FancyArrowPatch((40, 28), (40, 48), arrowstyle="<|-|>", mutation_scale=10,
                                 color="#333", lw=1.2, zorder=6))
    ax.text(41.5, 38, "height", ha="left", va="center", fontsize=10.5, rotation=90)
    arrow(ax, (50, 40), (57, 40), lw=1.8)
    seed_bag(ax, 64, 34, 8, 9, 14)
    ax.text(64, 30.5, "seeds", ha="center", va="top", fontsize=11.5)
    stopwatch(ax, 85, 42)
    ax.text(85, 33.5, "fast &\nlow cost", ha="center", va="top", fontsize=11, color=INTACT,
            fontweight="bold", linespacing=1.1)
    for k, (yr, mk) in enumerate((("2025", "o"), ("2026", "s"))):
        ax.text(70 + 15 * k, 15, yr, ha="center", va="center", fontsize=11,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#888"))


def panel_curve(ax):
    frame(ax, 2, "Seeds expected for size", "all plants → seed expected for plant size")
    x0, y0, w, h = 16, 18, 74, 66
    ax.plot([x0, x0, x0 + w], [y0 + h, y0, y0], color="#333", lw=1.5, zorder=3)
    xx = np.array([x0, x0 + w])
    line = lambda x, off=0: y0 + 12 + (x - x0) * 0.6 + off
    ax.fill_between(xx, line(xx, -9), line(xx, 9), color=INTACT, alpha=0.16, lw=0, zorder=1)
    ax.plot(xx, line(xx), color=INTACT, lw=2.4, zorder=3)
    ax.plot(xx, line(xx, -9), color=MATE_LIM, lw=1.6, ls="--", zorder=3)
    ax.plot(xx, line(xx, 9), color=OVER, lw=1.6, ls="--", zorder=3)
    rng = np.random.default_rng(3)
    px = rng.uniform(x0 + 3, x0 + w - 3, 70)
    ax.scatter(px, np.clip(line(px, rng.normal(0, 5.5, 70)), y0 + 1, y0 + h), s=14, c=GREY_PT,
               ec="white", lw=0.3, zorder=2)
    ax.text(x0 + w / 2, y0 - 3, "plant size (crown)", ha="center", va="top", fontsize=12)
    ax.text(x0 - 3, y0 + h / 2, "seeds per plant", ha="right", va="center", fontsize=12,
            rotation=90)
    ax.text(x0 + w - 1, line(x0 + w, 9) + 2, "twice", ha="right", va="bottom", fontsize=11,
            color=OVER)
    ax.text(x0 + w + 1, line(x0 + w) - 2, "expected", ha="right", va="top", fontsize=11,
            color=INTACT, fontweight="bold")
    ax.text(x0 + w - 1, line(x0 + w, -9) - 3, "half", ha="right", va="top", fontsize=11,
            color=MATE_LIM)


def panel_compare(ax):
    frame(ax, 3, "Compare each population", "two seasons: below, as expected or above")
    xs = {0.5: 32, 1: 55, 2: 78}                         # log-scaled positions
    for v, x in xs.items():
        ax.plot([x, x], [24, 84], color={0.5: MATE_LIM, 1: "#777", 2: OVER}[v], lw=1.4,
                ls="-" if v == 1 else "--", zorder=1)
        ax.text(x, 21, {0.5: "half", 1: "as expected", 2: "twice"}[v], ha="center", va="top",
                fontsize=11, color={0.5: MATE_LIM, 1: "#555", 2: OVER}[v])
    ax.add_patch(Rectangle((32, 24), 46, 60, fc=INTACT, alpha=0.10, lw=0, zorder=0))
    rows = [("pop. A", (57, 52), INTACT), ("pop. B", (21, 26), MATE_LIM), ("pop. C", (86, 83), OVER),
            ("pop. D", (50, 61), INTACT)]
    for k, (lab, (x25, x26), c) in enumerate(rows):
        y = 76 - 15 * k
        ax.plot([x25, x26], [y, y], color=c, lw=1.2, zorder=2)
        ax.scatter([x25], [y], s=110, marker="o", c=c, ec="white", lw=0.8, zorder=3)
        ax.scatter([x26], [y], s=110, marker="s", c=c, ec="white", lw=0.8, zorder=3)
        ax.text(13, y, lab, ha="right", va="center", fontsize=11.5)
    ax.scatter([18], [91], s=60, marker="o", c="#555"), ax.text(21, 91, "2025", va="center",
                                                                fontsize=11)
    ax.scatter([38], [91], s=60, marker="s", c="#555"), ax.text(41, 91, "2026", va="center",
                                                                fontsize=11)
    ax.text(55, 13.5, "seed set ÷ seed expected for size", ha="center", va="top", fontsize=11,
            color="#333")


def panel_candidates(ax):
    frame(ax, 4, "Candidate sources", "intact = candidate source → seed bulking")
    rows = [(76, INTACT, "intact", "source"), (49, MATE_LIM, "mate-limited", "target"),
            (22, OVER, "over-producing", "excluded")]
    for y, c, lab, role in rows:
        ax.add_patch(Ellipse((22, y), 28, 15, fc=SPOT, ec=c, lw=3, zorder=1))
        for dx in (-6, 6):
            lepa(ax, 22 + dx, y - 5, h=8, faded=c == MATE_LIM and dx > 0, fruits=c != MATE_LIM)
        ax.text(40, y + 2, lab, ha="left", va="bottom", fontsize=12.5, color=c, fontweight="bold")
        ax.text(40, y - 1, f"= {role}", ha="left", va="top", fontsize=12.5, color=c)
    cross_mark(ax, 22, 22, s=6, col=OVER, lw=2.5)
    ax.text(40, 64.5, "P34 · P35 · P43", ha="left", va="top", fontsize=10.5, color=INTACT,
            style="italic")
    arrow(ax, (62, 78), (70, 78), col=INTACT, lw=2.2)    # → seed bulking
    greenhouse(ax, 72, 70, 22, 18)
    for k in range(3):
        ax.add_patch(Rectangle((75.5 + 6 * k, 71), 3, 1.5, fc=rc.POT, ec=rc.POT_EDGE, lw=0.6,
                               zorder=2))
        lepa(ax, 77 + 6 * k, 72.5, h=4.5)
    ax.text(83, 66.5, "seed bulking", ha="center", va="top", fontsize=11, color=INTACT,
            fontweight="bold")
    dna(ax, 83, 40, h=11, w=2.4)
    ax.text(83, 37, "next: SRK\ncheck", ha="center", va="top", fontsize=10.5, color="#555",
            linespacing=1.1)
    arrow(ax, (83, 63), (83, 53), col="#888", lw=1.5, ls="--")


PANELS = (panel_measure, panel_curve, panel_compare, panel_candidates)


def draw(n_shown: int):
    rc.RNG = np.random.default_rng(11)                   # same plants in every frame
    fig, axes = plt.subplots(1, len(PANELS), figsize=(24, 7.0))
    for k, (ax, f) in enumerate(zip(axes, PANELS)):
        if k < n_shown:
            f(ax)
        else:
            ax.set_xlim(0, 100), ax.set_ylim(0, 100), ax.set_aspect("equal"), ax.axis("off")
            ax.set_title(" ", fontsize=18, pad=10)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.04, wspace=0.08)
    fig.canvas.draw()
    for k, (a, b) in enumerate(zip(axes[:-1], axes[1:])):
        if k + 1 >= n_shown:
            break
        p = fig.transFigure.inverted().transform(a.transData.transform((100, 50)))
        q = fig.transFigure.inverted().transform(b.transData.transform((0, 50)))
        fig.add_artist(FancyArrowPatch(p, q, transform=fig.transFigure, arrowstyle="-|>",
                                       mutation_scale=26, color="#333", lw=3))
    return fig


def main():
    OUT.mkdir(exist_ok=True)
    fig = draw(len(PANELS))
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"source_screening_cartoon.{ext}", dpi=200 if ext == "png" else None)
    plt.close(fig)
    for n in range(1, len(PANELS)):
        fig = draw(n)
        fig.savefig(OUT / f"source_screening_cartoon_step{n}.png", dpi=200)
        plt.close(fig)
    print(f"  wrote {OUT}/source_screening_cartoon.png (+ _step1–3 frames for slides)")


if __name__ == "__main__":
    main()
