"""restoration_program_cartoon.py — the population restoration program, with
slickspot peppergrass (Lepidium papilliferum) as the model (Seed_bulking
Figure 2). Replaces restoration_cartoon.py's figure, which is kept.

Objectives are told apart from techniques: each panel's title is the
**objective** (what), and a gray tag at its foot names the **technique**
(how).

  The challenge (very little seed to start with) is drawn here
     (panel_challenge) but shown at the top of Figure 1
     (breeding_system_cartoon.py), not in this figure.
  1. Prioritize populations for population restoration — reproductively
     intact populations are the sources (top conservation priority: losing
     one compromises the whole program) and go on to step 2; mate-limited
     populations are rescued (objective 3); populations under inbreeding
     depression (over-producing) are excluded as sources. Technique: seed
     set for plant size, drawn in Figure 3 (source_screening_cartoon.py).
  2. Multiply seed — technique: seed bulking in the greenhouse.
  3. Rescue mate-limited populations (header over two steps):
     3a. Breed seed with new S-alleles (limited genetic diversity in the
     target) — technique: informed breeding. Seedlings of the source
     (pollen donors) and the target (mothers) are grown together; source
     pollen goes on target flowers, pairs chosen by SRK so that they share
     no S-allele; the seed is harvested from the target mothers.
     3b. Sow the seed into the population — technique: seed outcasting;
     seed set and S-alleles monitored.

Icons and colors as restoration_cartoon.py; invented plants and alleles.
Output: figures_seed_bulking/restoration_program_cartoon.{png,pdf}, plus
restoration_program_cartoon_step1–3.png (same size, later parts blank) to
build a slide; landscape slides restoration_program_cartoon_slide1 (steps 1
and 2) and _slide2 (3a and 3b), {png,pdf}, with _step frames.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import (Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon,
                                Rectangle)
from matplotlib.path import Path as MPath

import restoration_cartoon as rc
from type_icons import add_icon, lifebuoy, no_sign
from restoration_cartoon import (ALLELE, INTACT, LEAF, MATE_LIM, NO, OUT, OVER, POT, POT_EDGE,
                                 SPOT, SPOT_EDGE, STEM, arrow, bee, brush, clipboard, cross_mark,
                                 dna, greenhouse, lepa, seed, seed_bag, tube)

W, H = 16.0, 17.85                                       # figure size, inches
TECH = "#546E7A"                                         # technique tags


def axes_in(fig, x, y_top, w, h):
    """Axes placed in inches, y measured from the top of the figure."""
    return fig.add_axes([x / W, 1 - (y_top + h) / H, w / W, h / H])


def frame(ax, title, n=None, xmax=100, ymax=100):
    ax.set_xlim(0, xmax), ax.set_ylim(0, ymax), ax.set_aspect("equal")
    ax.set_xticks([]), ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor("#BBBBBB")
    ax.set_title(f"{n}. {title}" if n else title, loc="left", fontsize=18, fontweight="bold",
                 pad=10)


def technique(ax, text, x=50, y=5):
    """Gray tag at the foot of a panel: the technique (how)."""
    ax.text(x, y, f"technique:  {text}", ha="center", va="center", fontsize=13.5, color="white",
            fontweight="bold", zorder=10,
            bbox=dict(boxstyle="round,pad=0.45", fc=TECH, ec="none"))


def shield(ax, x, y, s=1.0, col=INTACT):
    """Conservation-priority shield with a star."""
    pts = [[x - 7 * s, y + 8 * s], [x, y + 10 * s], [x + 7 * s, y + 8 * s], [x + 6 * s, y - 2 * s],
           [x, y - 9 * s], [x - 6 * s, y - 2 * s]]
    ax.add_patch(Polygon(pts, closed=True, fc=col, ec="#004D40", lw=1.6, zorder=6,
                         joinstyle="round"))
    ang = np.deg2rad(90 + 36 * np.arange(10))
    r = np.where(np.arange(10) % 2 == 0, 4.4 * s, 1.8 * s)
    ax.add_patch(Polygon(np.c_[x + r * np.cos(ang), y + 0.5 * s + r * np.sin(ang)], fc="white",
                         ec="none", zorder=7))


def magnifier(ax, x, y, r=8):
    ax.add_patch(Circle((x, y), r, fc="none", ec="#37474F", lw=2.4, zorder=8))
    a = np.deg2rad(-45)
    ax.plot([x + r * np.cos(a), x + 1.9 * r * np.cos(a)], [y + r * np.sin(a), y + 1.9 * r * np.sin(a)],
            color="#37474F", lw=5, solid_capstyle="round", zorder=8)


# ---- the challenge ---------------------------------------------------------------
def threatened_tag(ax, x, y, size=12.5, ha="left"):
    """Status tag: the species is federally listed as threatened."""
    ax.text(x, y, "federally threatened (ESA)", ha=ha, va="center", fontsize=size, color=NO,
            fontweight="bold", zorder=10,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=NO, lw=1.6))


def consequences(ax, x, y, dy, size, indent=None):
    """Why seed scarcity impedes restoration: two ✗ lines (the problems to solve)."""
    for k, line in enumerate(("far too little to restore\npopulations directly",
                              "collecting more would deplete\nthe populations it comes from")):
        ax.text(x, y - k * dy, "✗", ha="left", va="center", fontsize=size + 6, color=NO,
                fontweight="bold")
        ax.text(x + (0.5 * dy if indent is None else indent), y - k * dy, line, ha="left",
                va="center", fontsize=size, color=NO, fontweight="bold", linespacing=1.15)


def panel_challenge(ax, title="The challenge: very little seed to start with"):
    """Wide panel (200 × 38 data units): the top of Figure 1 (breeding_system_cartoon)."""
    frame(ax, title, xmax=200, ymax=38)
    threatened_tag(ax, 3, 34, size=13.5)
    ax.text(50, 34, "wild populations: limited fertility — few fruits", ha="left", va="center",
            fontsize=13.5, color="#555", style="italic")
    # slick spots jittered up and down, so they can be larger without touching
    spots = [(16, 12.5, 30, 20, 3), (39.5, 22, 22, 15.5, 2), (62.5, 12, 30, 20, 3),
             (89, 21, 22, 16.5, 2)]
    for x, y, w, h, n in spots:
        ax.add_patch(Ellipse((x, y), w, h, fc=SPOT, ec=SPOT_EDGE, lw=1.4, zorder=1))
        for dx in np.linspace(-(n - 1) * 4, (n - 1) * 4, n):
            lepa(ax, x + dx, y - 6, h=10, faded=dx != 0, fruits=1 if dx == 0 else 0)
    arrow(ax, (101.5, 16), (108.5, 16), lw=2.6)
    seed_bag(ax, 115, 9, 9, 10, 3)
    ax.text(122, 16, "very little\nseed", ha="left", va="center", fontsize=16, fontweight="bold",
            linespacing=1.15)
    consequences(ax, 141, 22, 12, 15)


def panel_challenge_slide(ax):
    """The challenge alone, filling a 16:9 slide (160 × 90 data units, no frame):
    the same content as panel_challenge, laid out for the talk."""
    ax.set_xlim(0, 160), ax.set_ylim(0, 90), ax.set_aspect("equal"), ax.axis("off")
    ax.text(4, 83, "The challenge: very little seed to start with — why?", ha="left", va="center",
            fontsize=34, fontweight="bold")
    threatened_tag(ax, 5, 72, size=22)
    ax.text(62, 72, "wild populations: limited fertility — few fruits", ha="left", va="center",
            fontsize=22, color="#555", style="italic")
    spots = [(26, 46, 44, 28, 3), (71, 51, 34, 22, 2), (29, 16, 34, 22, 2), (71, 18, 44, 28, 3)]
    for x, y, w, h, n in spots:
        ax.add_patch(Ellipse((x, y), w, h, fc=SPOT, ec=SPOT_EDGE, lw=1.8, zorder=1))
        for dx in np.linspace(-(n - 1) * 6, (n - 1) * 6, n):
            lepa(ax, x + dx, y - 7, h=15, faded=dx != 0, fruits=1 if dx == 0 else 0)
    arrow(ax, (92, 50), (102, 50), lw=4)
    seed_bag(ax, 113, 43, 13, 14, 3)
    ax.text(124, 50, "very little\nseed", ha="left", va="center", fontsize=30, fontweight="bold",
            linespacing=1.15)
    consequences(ax, 98, 28, 17, 21, indent=6)


# ---- 1. identify source populations ------------------------------------------------
def badge(ax, kind, x, y, size):
    """Type badge of type_icons.py centered at (x, y), `size` data units wide
    (panel axes: square, data 0–100)."""
    pos = ax.get_position()
    w = size / 100 * pos.width
    h = size / 100 * pos.height
    cx, cy = pos.x0 + x / 100 * pos.width, pos.y0 + y / 100 * pos.height
    add_icon(ax.figure, kind, [cx - w / 2, cy - h / 2, w, h])


def panel_sources(ax):
    """Prioritize populations: intact = source (priority 1, protect) → step 2;
    mate-limited → rescue (objective 3); over-producing (inbreeding
    depression) → excluded. Badges as in Figure 3; action symbols: shield,
    ring buoy (rescue), no-sign (avoid)."""
    frame(ax, "Prioritize populations for population restoration", 1)
    # reproductively intact = source: top conservation priority, on to step 2
    badge(ax, "Reproductively intact", 58, 70, 24)
    magnifier(ax, 58, 70, r=14.5)
    ax.text(58, 89.5, "reproductively intact (source)", ha="center", fontsize=13.5, color=INTACT,
            fontweight="bold")
    shield(ax, 17, 72, s=1.05)
    ax.text(17, 56.5, "top conservation\npriority: protect", ha="center", va="top",
            fontsize=11.5, color=INTACT, fontweight="bold", linespacing=1.1)
    ax.plot([4, 96], [47, 47], color="#CCCCCC", lw=1.2, ls=(0, (4, 3)), zorder=0)
    # the other two types
    for y, col, kind, head, action in (
            (35, MATE_LIM, "Mate-limited", "mate-limited", "→ rescue (objective 3)"),
            (16, OVER, "Over-producing", "over-producing",
             "→ inbreeding depression:\n    excluded as source")):
        badge(ax, kind, 17, y, 17)
        ax.text(31, y + 2.5, head, ha="left", va="bottom", fontsize=13.5, color=col,
                fontweight="bold")
        ax.text(31, y - 0.5, action, ha="left", va="top", fontsize=12.5, color=col,
                linespacing=1.1)
    lifebuoy(ax, 71, 38, r=5)                          # mate-limited: rescue
    ax.text(71, 31.5, "rescue", ha="center", va="top", fontsize=11.5, color=MATE_LIM,
            fontweight="bold")
    no_sign(ax, 71, 20, r=5)                              # over-producing: avoid as source
    ax.text(71, 13.5, "avoid", ha="center", va="top", fontsize=11.5, color="#C62828",
            fontweight="bold")
    technique(ax, "seed set for plant size (see Figure 3)", y=4)


def panel_sources_talk(ax, use=True):
    """Talk version of panel_sources, no title: the three population types
    aligned, one row each (badge, type); use=True adds what each means for
    restoration (use of the type, action symbol). Rows do not move between
    the two versions."""
    frame(ax, "")
    rows = ((85, INTACT, "Reproductively intact", "reproductively intact",
             "→ source of seed: top\n    conservation priority"),
            (55.5, MATE_LIM, "Mate-limited", "mate-limited", "→ genetic rescue"),
            (26, OVER, "Over-producing", "over-producing",
             "→ inbreeding depression:\n    excluded as source"))
    for y, col, kind, head, action in rows:              # rows fill the panel above the tag
        badge(ax, kind, 15, y, 24)
        ax.text(31, y + 2.5, head, ha="left", va="bottom", fontsize=19, color=col,
                fontweight="bold")
        if use:
            ax.text(31, y - 0.5, action, ha="left", va="top", fontsize=17, color=col,
                    linespacing=1.1)
    for y in (70.25, 40.75):
        ax.plot([4, 96], [y, y], color="#CCCCCC", lw=1.2, ls=(0, (4, 3)), zorder=0)
    technique(ax, "seed set for plant size", y=5)
    if not use:
        return
    shield(ax, 86, 88, s=0.72)
    ax.text(86, 76.5, "protect", ha="center", va="center", fontsize=15, color=INTACT,
            fontweight="bold")
    lifebuoy(ax, 86, 58.5, r=6.5)
    ax.text(86, 47, "rescue", ha="center", va="center", fontsize=15, color=MATE_LIM,
            fontweight="bold")
    no_sign(ax, 86, 29, r=6.5)
    ax.text(86, 17.5, "avoid", ha="center", va="center", fontsize=15, color="#C62828",
            fontweight="bold")


# ---- 2. multiply seed -----------------------------------------------------------
def panel_bulking(ax):
    frame(ax, "Multiply seed", 2)
    seed_bag(ax, 12, 66, 8, 8, 4, fc="#E0F2F1")
    ax.text(12, 62, "from\nsources", ha="center", va="top", fontsize=11.5, color=INTACT,
            linespacing=1.05)
    arrow(ax, (18, 70), (26, 70))
    greenhouse(ax, 28, 40, 66, 52)
    for r, yy in enumerate((65, 53, 41)):
        ax.text(31, yy + 3.5, f"M{r + 1}", fontsize=11, color="#555", va="center")
        for k in range(5):
            xx = 42 + 11.5 * k
            ax.add_patch(Rectangle((xx - 3, yy), 6, 2.5, fc=POT, ec=POT_EDGE, lw=0.8, zorder=2))
            lepa(ax, xx, yy + 2.5, h=7.5)
    brush(ax, 70, 77, s=0.9)
    arrow(ax, (61, 38), (61, 31))
    for x in (34, 45, 56, 67, 78, 89):
        seed_bag(ax, x, 16, 7, 8, 12)
    ax.text(14, 21, "many\nseeds", ha="center", va="center", fontsize=13, fontweight="bold",
            linespacing=1.1)
    technique(ax, "seed bulking (greenhouse)")


# ---- 3. genetic rescue ---------------------------------------------------------
SOURCE = [(0, 3), (2, 4), (5, 6)]                        # many S-alleles; plant 1 shares S1
TARGET = [(0, 1), (0, 1), (1, 0)]                        # few S-alleles


def panel_rescue(ax):
    """Seedlings of source and target grown together; source pollen on target
    flowers, pairs chosen by SRK (no shared S-allele); seed from the targets."""
    frame(ax, "Breed seed with new S-alleles", "3a")
    ax.text(50, 95, "source pollen → target flowers, pairs chosen by SRK", ha="center",
            fontsize=12.5, color="#333", style="italic")
    ax.add_patch(Ellipse((8, 84), 5, 2.6, angle=30, fc=LEAF, ec=STEM, lw=0.8))  # leaf → tube → DNA
    arrow(ax, (12, 84), (17, 84), lw=1.5)
    tube(ax, 21, 79, s=0.75)
    arrow(ax, (25, 84), (30, 84), lw=1.5)
    dna(ax, 35, 78, h=11, w=2.4)
    for k in range(8):                                   # S-allele key
        ax.add_patch(Circle((47 + 6.3 * k, 87), 2.1, fc=ALLELE[k], ec="k", lw=0.6))
        ax.text(47 + 6.3 * k, 83.5, f"S{k + 1}", ha="center", va="top", fontsize=9)
    greenhouse(ax, 22, 15, 74, 60)                       # seedlings of both, grown together
    xs = (38, 60, 82)
    for y, col, lab, sub, alls in ((47, INTACT, "source", "(pollen\ndonors)", SOURCE),
                                   (23, MATE_LIM, "target", "(mothers)", TARGET)):
        ax.text(12, y + 6, lab, ha="center", va="center", fontsize=13.5, color=col,
                fontweight="bold")
        ax.text(12, y + 1.5, sub, ha="center", va="top", fontsize=11, color=col,
                linespacing=1.05)
        for x, al in zip(xs, alls):
            ax.add_patch(Rectangle((x - 5, y - 4.2), 10, 1.2, fc=col, ec="none", zorder=2))
            lepa(ax, x, y, h=10, alleles=list(al), fruits=col == INTACT)
    for i, j in ((1, 0), (2, 1), (2, 2)):                # compatible: no shared S-allele
        ax.add_patch(FancyArrowPatch((xs[i] - 2, 55), (xs[j] + 2, 33), arrowstyle="-|>",
                                     mutation_scale=14, color=INTACT, lw=2.2,
                                     connectionstyle="arc3,rad=0.15", zorder=6))
    ax.add_patch(FancyArrowPatch((xs[0] - 4, 56), (xs[0] - 4, 33), arrowstyle="-|>",
                                 mutation_scale=14, color=NO, lw=2, ls=(0, (3, 2)),
                                 connectionstyle="arc3,rad=0.55", zorder=6))
    cross_mark(ax, xs[0] - 10.5, 44.5, s=1.8)           # shares S1 with the target: no cross
    brush(ax, 71, 42, s=0.8, angle=-35)
    arrow(ax, (84, 17), (88, 13))
    seed_bag(ax, 89, 4.5, 8, 7, 16, seed_cols=[ALLELE[k] for k in (0, 1, 2, 4, 5, 6)])
    ax.text(83, 8.5, "seed with\nnew S-alleles", ha="right", va="center", fontsize=11,
            linespacing=1.05)
    technique(ax, "informed breeding (SRK)", x=33, y=4)


# ---- 4. restore populations ------------------------------------------------------
def panel_restore(ax):
    frame(ax, "Sow the seed into the population", "3b")
    for x0, ok, al_sets in ((25, False, [(0, 1), (0, 1), (1, 0)]),
                            (75, True, [(0, 1), (2, 3), (4, 5), (6, 7), (0, 6)])):
        col = INTACT if ok else MATE_LIM
        ax.add_patch(Ellipse((x0, 48), 42, 32, fc=SPOT, ec=col, lw=3, zorder=1))
        ax.text(x0, 67, "after" if ok else "before: mate-limited", ha="center", fontsize=14,
                fontweight="bold", color=col)
        xs = np.linspace(x0 - 13, x0 + 13, len(al_sets))
        for k, (xx, al) in enumerate(zip(xs, al_sets)):
            lepa(ax, xx, 46 + (4 if k % 2 else -4), h=8, alleles=list(al), fruits=ok)
        if ok:                                           # bee flights between compatible plants
            for a, b in ((0, 1), (1, 2), (2, 3), (3, 4)):
                arrow(ax, (xs[a], 58), (xs[b], 58), col="#555", lw=1.2, rad=-0.35, ls="--",
                      head="-")
            bee(ax, x0 + 3, 61)
            seed_bag(ax, x0 + 14, 20, 7, 7, 14)
        else:
            bee(ax, x0, 61)
            seed_bag(ax, x0 + 12, 20, 7, 7, 2)
    arrow(ax, (46.5, 48), (53.5, 48), lw=2.5)
    seed_bag(ax, 50, 80, 7, 7, 12, seed_cols=ALLELE, fc="#E0F2F1")       # bred seed …
    for k, (sx, sy) in enumerate(((44, 75), (47, 71), (51, 73), (54, 69), (57, 74), (60, 70),
                                  (63, 66), (66, 71))):                   # … cast onto the site
        seed(ax, sx + 8, sy - 4, 1.1, 40 * k, ALLELE[k % len(ALLELE)])
    arrow(ax, (57, 81), (70, 70), col=INTACT, rad=-0.3, lw=1.8)
    clipboard(ax, 10, 15, s=1.0)
    ax.text(10, 12, "monitor", ha="center", va="top", fontsize=11, color="#555")
    technique(ax, "seed outcasting", x=55)


PANELS = (panel_sources, panel_bulking, panel_rescue, panel_restore)


def rescue_header(fig, ax_a, ax_b):
    """Objective 3 over its two steps (3a, 3b): title and a bracket."""
    hin = fig.get_figheight()
    to_fig = lambda ax, xy: fig.transFigure.inverted().transform(ax.transData.transform(xy))
    (xl, yt), (xr, _) = to_fig(ax_a, (0, 100)), to_fig(ax_b, (100, 100))
    yb = yt + 0.62 / hin                                 # bracket, above the step titles
    tick = 0.12 / hin
    fig.add_artist(plt.Line2D([xl, xl, xr, xr], [yb - tick, yb, yb, yb - tick],
                              transform=fig.transFigure, color="#333", lw=2.2))
    fig.text(xl, yb + 0.1 / hin, "3. Rescue mate-limited populations", ha="left", va="bottom",
             fontsize=18, fontweight="bold")              # same style as titles 1 and 2


def draw(n_shown: int):
    """n_shown: 1–4 = steps shown (4 = complete figure)."""
    rc.RNG = np.random.default_rng(7)                    # same plants in every frame
    fig = plt.figure(figsize=(W, H))
    side, gap = 7.2, 0.8
    xs, ys = (0.4, 0.4 + side + gap), (0.75, 0.75 + side + 2.4)
    grid = [[axes_in(fig, x, y, side, side) for x in xs] for y in ys]
    for k, (ax, f) in enumerate(zip([a for row in grid for a in row], PANELS)):
        if k < n_shown:
            f(ax)
        else:
            ax.axis("off")
    # connectors: 1 → 2, 2 → 3 (elbow), 3 → 4
    fig.canvas.draw()
    to_fig = lambda ax, xy: fig.transFigure.inverted().transform(ax.transData.transform(xy))
    style = dict(transform=fig.transFigure, arrowstyle="-|>", mutation_scale=26, color="#333",
                 lw=3)
    if n_shown >= 2:                                     # 1 → 2 leaves from the source
        fig.add_artist(FancyArrowPatch(to_fig(grid[0][0], (75, 70)), to_fig(grid[0][1], (0, 70)),
                                       shrinkA=4, shrinkB=6, **{**style, "color": INTACT}))
    if n_shown >= 3:
        rescue_header(fig, grid[1][0], grid[1][1] if n_shown >= 4 else grid[1][0])
    if n_shown >= 4:
        fig.add_artist(FancyArrowPatch(to_fig(grid[1][0], (100, 50)), to_fig(grid[1][1], (0, 50)),
                                       shrinkA=6, shrinkB=6, **style))
    if n_shown >= 3:
        (x0, y0), (x1, y1) = to_fig(grid[0][1], (50, 0)), to_fig(grid[1][0], (88, 100))
        y_mid = y0 - 0.012                               # ends above the objective-3 header
        y_end = y1 + 1.2 / H
        fig.add_artist(FancyArrowPatch(path=MPath([(x0, y0 - 0.002), (x0, y_mid), (x1, y_mid),
                                                  (x1, y_end)]), **style))
    return fig


def draw_slide(panels, n_shown: int):
    """Landscape slide: square panels side by side, arrows between them."""
    rc.RNG = np.random.default_rng(7)
    side, gap, m = 7.2, 1.0, 0.35
    extra = 0.75 if panel_rescue in panels else 0.0      # room for the objective-3 header
    w, h = 2 * m + len(panels) * side + (len(panels) - 1) * gap, 2 * m + 0.6 + side + extra
    fig = plt.figure(figsize=(w, h))
    axes = [fig.add_axes([(m + k * (side + gap)) / w, m / h, side / w, side / h])
            for k in range(len(panels))]
    for k, (ax, f) in enumerate(zip(axes, panels)):
        if k < n_shown:
            f(ax)
        else:
            ax.axis("off")
    fig.canvas.draw()
    to_fig = lambda ax, xy: fig.transFigure.inverted().transform(ax.transData.transform(xy))
    style = dict(transform=fig.transFigure, arrowstyle="-|>", mutation_scale=26, color="#333",
                 lw=3)
    if panel_rescue in panels and n_shown > panels.index(panel_rescue):
        i = panels.index(panel_rescue)
        rescue_header(fig, axes[i], axes[i + 1] if n_shown > i + 1 else axes[i])
    for k in range(min(n_shown, len(panels)) - 1):
        yy = 70 if panels[k + 1] is panel_bulking else 50
        a, b = axes[k], axes[k + 1]
        fig.add_artist(FancyArrowPatch(to_fig(a, (100, yy)), to_fig(b, (0, yy)), shrinkA=6,
                                       shrinkB=6, **style))
    return fig


SLIDES = {"slide1": (panel_sources, panel_bulking),     # objectives 1 – 2: sources, more seed
          "slide2": (panel_rescue, panel_restore)}     # objective 3: rescue (3a, 3b)


def draw_types(use: bool):
    """Talk: the three types alone (panel_sources_talk), one square panel."""
    rc.RNG = np.random.default_rng(7)
    side, m = 7.2, 0.1
    fig = plt.figure(figsize=(side + 2 * m, side + 2 * m))
    k = m / (side + 2 * m)
    panel_sources_talk(fig.add_axes([k, k, 1 - 2 * k, 1 - 2 * k]), use)
    return fig


def main():
    OUT.mkdir(exist_ok=True)
    for name, use in (("types", False), ("types_use", True)):
        fig = draw_types(use)
        fig.savefig(OUT / f"restoration_program_cartoon_{name}.png", dpi=200)
        plt.close(fig)
    fig = draw(4)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"restoration_program_cartoon.{ext}", dpi=200 if ext == "png" else None)
    plt.close(fig)
    for n in range(1, 4):                                # talk build, step by step
        fig = draw(n)
        fig.savefig(OUT / f"restoration_program_cartoon_step{n}.png", dpi=200)
        plt.close(fig)
    for name, panels in SLIDES.items():                  # landscape slides, built step by step
        fig = draw_slide(panels, len(panels))
        for ext in ("png", "pdf"):
            fig.savefig(OUT / f"restoration_program_cartoon_{name}.{ext}",
                        dpi=200 if ext == "png" else None)
        plt.close(fig)
        for n in range(1, len(panels)):
            fig = draw_slide(panels, n)
            fig.savefig(OUT / f"restoration_program_cartoon_{name}_step{n}.png", dpi=200)
            plt.close(fig)
    print(f"  wrote {OUT}/restoration_program_cartoon.png (+ _step1–3 frames) and the "
          f"landscape slides _slide1 (1, 2) and _slide2 (3a, 3b)")


if __name__ == "__main__":
    main()
