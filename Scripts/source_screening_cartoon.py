"""source_screening_cartoon.py — our approach to sort populations for
population restoration, with slickspot peppergrass (Lepidium papilliferum)
as the model (Seed_bulking Figure 3). It is the technique of step 1 of
restoration_program_cartoon.py (Figure 3).

Four steps, left to right:
  1. **Plant phenotyping** — plant size (crown width, height) and seed
     yield of every flowering plant, in two seasons: fast and low cost, no
     genotyping.
  2. **Predicting reproductive status** — (a) build: all plants together give the seed yield a
     plant of a given size is expected to set (band: half to twice);
     (b) deploy by population: each population's seed yield relative to
     that expectation, per season.
  3. **Sort populations for population restoration** — reproductively
     intact = source (seed bulking; top conservation priority);
     mate-limited = genetic rescue (flaw: reduced fecundity);
     over-producing = excluded as source (flaw: inbreeding depression).
     The sort is a set of hypotheses, tested in panel 4.
  4. **Test with genetics** — why (resource or pollinator availability can
     also drive seed yield; our research shows pollinators are not
     limiting), and what the SRK gene (S-types) should show for each type.

Same icons and colors as restoration_cartoon.py; invented plants and
values. Output: figures_seed_bulking/source_screening_cartoon.{png,pdf}
(landscape, for the talk) plus source_screening_cartoon_step1–3.png (same
size, later steps blank) to build a slide, and
source_screening_cartoon_compact.{png,pdf} (2 × 2 grid: the web document,
Seed_bulking Figure 4).
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path as MPath
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle

import restoration_cartoon as rc
from type_icons import add_icon, lifebuoy, no_sign
from restoration_cartoon import (INTACT, MATE_LIM, OUT, OVER, SPOT, SPOT_EDGE, arrow, cross_mark,
                                 dna, lepa, seed_bag)

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


def shield(ax, x, y, s=1.0, col=INTACT):
    pts = [[x - 7 * s, y + 8 * s], [x, y + 10 * s], [x + 7 * s, y + 8 * s], [x + 6 * s, y - 2 * s],
           [x, y - 9 * s], [x - 6 * s, y - 2 * s]]
    ax.add_patch(Polygon(pts, closed=True, fc=col, ec="#004D40", lw=1.4, zorder=6))
    ang = np.deg2rad(90 + 36 * np.arange(10))
    r = np.where(np.arange(10) % 2 == 0, 4.4 * s, 1.8 * s)
    ax.add_patch(Polygon(np.c_[x + r * np.cos(ang), y + 0.5 * s + r * np.sin(ang)], fc="white",
                         ec="none", zorder=7))


def badge(ax, kind, x, y, size):
    """Type badge of type_icons.py centered at (x, y), `size` data units wide
    (panel axes: square, data 0–100)."""
    pos = ax.get_position()
    cx, cy = pos.x0 + x / 100 * pos.width, pos.y0 + y / 100 * pos.height
    w, h = size / 100 * pos.width, size / 100 * pos.height
    add_icon(ax.figure, kind, [cx - w / 2, cy - h / 2, w, h])


KINDS = {INTACT: "Reproductively intact", MATE_LIM: "Mate-limited", OVER: "Over-producing"}


# ---- 1. plant phenotyping ------------------------------------------------------------
def panel_phenotyping(ax):
    rc.frame(ax, 1, "Plant phenotyping", "plant size + seed yield, every plant, two seasons")
    ax.add_patch(Ellipse((50, 74), 70, 22, fc=SPOT, ec=SPOT_EDGE, lw=1.4, zorder=1))
    for dx, h in ((-26, 7), (-14, 10), (-2, 6), (10, 12), (22, 8)):
        lepa(ax, 50 + dx, 67, h=h)
    ax.text(80, 62.5, "slick spot", ha="center", va="top", fontsize=11.5, color="#8D6E3F",
            fontweight="bold", style="italic")
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
    ax.text(25, 20.5, "plant size", ha="center", va="top", fontsize=12, fontweight="bold")
    arrow(ax, (50, 40), (57, 40), lw=1.8)
    seed_bag(ax, 64, 34, 8, 9, 14)
    ax.text(64, 30.5, "seed yield", ha="center", va="top", fontsize=12, fontweight="bold")
    stopwatch(ax, 85, 42)
    ax.text(85, 33.5, "fast &\nlow cost", ha="center", va="top", fontsize=11, color=INTACT,
            fontweight="bold", linespacing=1.1)
    for k, yr in enumerate(("2025", "2026")):
        ax.text(70 + 15 * k, 15, yr, ha="center", va="center", fontsize=11,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#888"))


# ---- 2. modeling: build, then deploy by population -----------------------------------------
def panel_modeling(ax):
    rc.frame(ax, 2, "Predicting reproductive status", "build on all plants → deploy by population")
    # (a) build: seed yield expected for plant size
    ax.text(3, 95, "a  build", ha="left", va="center", fontsize=13, fontweight="bold")
    ax.text(25, 95, "all plants → seed yield expected for size", ha="left", va="center",
            fontsize=11, color="#555")
    x0, y0, w, h = 18, 58, 74, 31
    ax.plot([x0, x0, x0 + w], [y0 + h, y0, y0], color="#333", lw=1.3, zorder=3)
    xx = np.array([x0, x0 + w])
    line = lambda x, off=0: y0 + 5 + (x - x0) * 0.3 + off
    ax.fill_between(xx, line(xx, -4.5), line(xx, 4.5), color=INTACT, alpha=0.16, lw=0, zorder=1)
    ax.plot(xx, line(xx), color=INTACT, lw=2.2, zorder=3)
    ax.plot(xx, line(xx, -4.5), color=MATE_LIM, lw=1.3, ls="--", zorder=3)
    ax.plot(xx, line(xx, 4.5), color=OVER, lw=1.3, ls="--", zorder=3)
    rng = np.random.default_rng(3)
    px = rng.uniform(x0 + 2, x0 + w - 2, 55)
    ax.scatter(px, np.clip(line(px, rng.normal(0, 2.8, 55)), y0 + 1, y0 + h), s=9, c=GREY_PT,
               ec="white", lw=0.3, zorder=2)
    ax.text(x0 + w / 2, y0 - 1.5, "plant size", ha="center", va="top", fontsize=10.5)
    ax.text(x0 - 2, y0 + h / 2, "seed yield", ha="right", va="center", fontsize=10.5,
            rotation=90)
    ax.text(x0 + w - 1, line(x0 + w - 1) - 2.5, "expected", ha="right", va="top", fontsize=10.5,
            color=INTACT, fontweight="bold", zorder=5,
            bbox=dict(boxstyle="round,pad=0.15", fc="#D9EEEA", ec="none"))
    # (b) deploy by population: seed yield ÷ expected, per season
    arrow(ax, (7, 88), (7, 49), lw=2)                 # build → deploy
    ax.text(3, 45, "b  deploy", ha="left", va="center", fontsize=13, fontweight="bold")
    ax.text(27, 45, "each population, each season", ha="left", va="center", fontsize=11,
            color="#555")
    xs = {0.5: 36, 1: 57, 2: 78}                         # log-scaled positions
    ytop, ybot = 39, 15
    ax.add_patch(Rectangle((36, ybot), 42, ytop - ybot, fc=INTACT, alpha=0.10, lw=0, zorder=0))
    for v, x in xs.items():
        col = {0.5: MATE_LIM, 1: "#777", 2: OVER}[v]
        ax.plot([x, x], [ybot, ytop], color=col, lw=1.2, ls="-" if v == 1 else "--", zorder=1)
        ax.text(x, ybot - 1.5, {0.5: "half", 1: "expected", 2: "twice"}[v], ha="center",
                va="top", fontsize=10, color=col)
    for k, (lab, (x25, x26), c) in enumerate((("pop. A", (59, 54), INTACT),
                                              ("pop. B", (24, 28), MATE_LIM),
                                              ("pop. C", (87, 84), OVER))):
        y = 34 - 8.5 * k
        ax.plot([x25, x26], [y, y], color=c, lw=1.1, zorder=2)
        ax.scatter([x25], [y], s=70, marker="o", c=c, ec="white", lw=0.7, zorder=3)
        ax.scatter([x26], [y], s=70, marker="s", c=c, ec="white", lw=0.7, zorder=3)
        ax.text(16, y, lab, ha="right", va="center", fontsize=10.5)
        badge(ax, KINDS[c], 4.5, y, 7.5)                  # type, as in panel 3
    ax.scatter([86], [44.5], s=35, marker="o", c="#555")
    ax.text(88.5, 44.5, "2025", va="center", fontsize=9.5)
    ax.scatter([86], [40.5], s=35, marker="s", c="#555")
    ax.text(88.5, 40.5, "2026", va="center", fontsize=9.5)


# ---- 3. sort populations for population restoration ----------------------------------------
def panel_sort(ax):
    """Sort populations: the three types and what to do with each — hypotheses
    from phenotyping + prediction, tested with genetics in panel 4."""
    rc.frame(ax, 3, "Sort populations for restoration", "")
    rows = [(83, INTACT, "reproductively intact", "→ source for seed bulking", None),
            (57, MATE_LIM, "mate-limited", "→ genetic rescue", "flaw: reduced fecundity"),
            (31, OVER, "over-producing", "→ excluded as source", "flaw: inbreeding depression")]
    for y, c, lab, role, flaw in rows:
        badge(ax, KINDS[c], 15, y, 19)                   # type badges, as in Figures 2 – 3
        ax.text(31, y + 3, lab, ha="left", va="bottom", fontsize=13, color=c, fontweight="bold")
        ax.text(31, y + 0.5, role, ha="left", va="top", fontsize=12.5, color=c)
        if flaw:
            ax.text(31, y - 5.5, flaw, ha="left", va="top", fontsize=11, color=c, style="italic")
    # what to do with each type: protect, rescue, avoid
    for (y, c, *_), lab in zip(rows, ("top priority", "rescue", "avoid")):
        ax.text(88, y - 8.5, lab, ha="center", va="top", fontsize=10.5, color=c,
                fontweight="bold")
    shield(ax, 88, 84, s=0.7)
    lifebuoy(ax, 88, 57, r=5.5)
    no_sign(ax, 88, 31, r=5)
    ax.text(50, 8, "hypotheses → tested with genetics (panel 4)", ha="center", va="center",
            fontsize=12.5, color="#222", style="italic",
            bbox=dict(boxstyle="round,pad=0.45", fc="#F3F4F6", ec="#888", lw=1.1))


# ---- 4. test with genetics ----------------------------------------------------------------
GEN_ROWS = (  # type color, what genetics should show if the call is right
    (INTACT, "reproductively intact", "many S-types, many compatible mates"),
    (MATE_LIM, "mate-limited", "few S-types, few compatible mates"),
    (OVER, "over-producing", "seeds carry only the mother's S-types (selfing)"),
)


def panel_genetics(ax):
    """Panel 4: why the sort must be tested (resource or pollinator availability
    can also drive seed yield; pollinators are not limiting), and what the
    SRK gene should show for each type."""
    rc.frame(ax, 4, "Test with genetics", "")
    ax.add_patch(FancyBboxPatch((3, 70), 94, 26, boxstyle="round,pad=0,rounding_size=2",
                                fc="#F3F4F6", ec="#888", lw=1.2, zorder=1))
    ax.text(7, 91, "Why test?", ha="left", va="center", fontsize=13, fontweight="bold",
            color="#222", zorder=2)
    ax.text(7, 85, "seed yield can also be driven by", ha="left", va="center", fontsize=12.5,
            color="#222", zorder=2)
    ax.text(7, 80, "resource or pollinator availability", ha="left", va="center",
            fontsize=12.5, color="#222", zorder=2)
    ax.text(7, 74.5, "our research: pollinators are not limiting → genetics", ha="left",
            va="center", fontsize=11.5, color="#222", style="italic", zorder=2)
    # leaf → tube → DNA: the SRK gene gives the S-types of every plant
    ax.add_patch(Ellipse((9, 60), 5, 2.6, angle=30, fc=rc.LEAF, ec=rc.STEM, lw=0.8))
    arrow(ax, (13, 60), (18, 60), lw=1.5)
    rc.tube(ax, 22, 55, s=0.75)
    arrow(ax, (26, 60), (31, 60), lw=1.5)
    dna(ax, 35, 54, h=11, w=2.4)
    ax.text(42, 62.5, "SRK gene:", ha="left", va="center", fontsize=12.5, fontweight="bold",
            color="#222")
    ax.text(42, 57.5, "the S-types of every plant", ha="left", va="center", fontsize=12.5,
            color="#222")
    ax.text(50, 48.5, "what genetics should show", ha="center", va="center", fontsize=11.5,
            color="#555", style="italic")
    for (c, lab, sign), y in zip(GEN_ROWS, (38, 24, 10)):
        badge(ax, KINDS[c], 11, y, 13)
        ax.text(21, y + 2.6, lab, ha="left", va="center", fontsize=12.5, color=c,
                fontweight="bold")
        ax.text(21, y - 2.6, sign, ha="left", va="center", fontsize=11, color="#222")


PANELS = (panel_phenotyping, panel_modeling, panel_sort, panel_genetics)


def draw(n_shown: int):
    """Four square panels in a row (landscape, slides), arrows between them."""
    rc.RNG = np.random.default_rng(11)                   # same plants in every frame
    side, gap, m = 7.2, 1.0, 0.35
    n = len(PANELS)
    w, h = 2 * m + n * side + (n - 1) * gap, 2 * m + 0.6 + side
    fig = plt.figure(figsize=(w, h))
    axes = [fig.add_axes([(m + k * (side + gap)) / w, m / h, side / w, side / h]) for k in range(n)]
    for k, (ax, f) in enumerate(zip(axes, PANELS)):
        if k < n_shown:
            f(ax)
        else:
            ax.axis("off")
    fig.canvas.draw()
    to_fig = lambda ax, xy: fig.transFigure.inverted().transform(ax.transData.transform(xy))
    for k in range(min(n_shown, n) - 1):
        fig.add_artist(FancyArrowPatch(to_fig(axes[k], (100, 50)), to_fig(axes[k + 1], (0, 50)),
                                       transform=fig.transFigure, arrowstyle="-|>",
                                       mutation_scale=26, color="#333", lw=3, shrinkA=6,
                                       shrinkB=6))
    return fig


def draw_compact():
    """The web document: the four panels as a 2 × 2 grid (1 → 2, then down to
    3 → 4), so the whole figure fits on a screen with larger text."""
    rc.RNG = np.random.default_rng(11)
    side, gap, m, row = 7.2, 1.0, 0.35, 1.3
    w, h = 2 * m + 2 * side + gap, 2 * m + 2 * (0.6 + side) + row - 0.6
    fig = plt.figure(figsize=(w, h))
    ys = (m + side + row, m)                              # bottom of the top and lower rows
    xs = (m, m + side + gap)
    axes = [fig.add_axes([x / w, y / h, side / w, side / h]) for y in ys for x in xs]
    for ax, f in zip(axes, PANELS):
        f(ax)
    fig.canvas.draw()
    to_fig = lambda ax, xy: fig.transFigure.inverted().transform(ax.transData.transform(xy))
    style = dict(transform=fig.transFigure, arrowstyle="-|>", mutation_scale=26, color="#333", lw=3)
    for a, b in ((0, 1), (2, 3)):
        fig.add_artist(FancyArrowPatch(to_fig(axes[a], (100, 50)), to_fig(axes[b], (0, 50)),
                                       shrinkA=6, shrinkB=6, **style))
    (x0, y0), (x1, y1) = to_fig(axes[1], (50, 0)), to_fig(axes[2], (88, 100))
    y_mid = y0 - 0.3 / h                                 # 2 → 3: elbow in the gap between rows
    fig.add_artist(FancyArrowPatch(path=MPath([(x0, y0 - 0.002), (x0, y_mid), (x1, y_mid),
                                              (x1, y1 + 0.08 / h)]), **style))
    return fig


def main():
    OUT.mkdir(exist_ok=True)
    for name, make in (("source_screening_cartoon", lambda: draw(len(PANELS))),
                       ("source_screening_cartoon_compact", draw_compact)):
        fig = make()
        for ext in ("png", "pdf"):
            fig.savefig(OUT / f"{name}.{ext}", dpi=200 if ext == "png" else None)
        plt.close(fig)
    for n in range(1, len(PANELS)):
        fig = draw(n)
        fig.savefig(OUT / f"source_screening_cartoon_step{n}.png", dpi=200)
        plt.close(fig)
    print(f"  wrote {OUT}/source_screening_cartoon.png (landscape, talk; + _step1–3 frames) "
          f"and _compact (web document)")


if __name__ == "__main__":
    main()
