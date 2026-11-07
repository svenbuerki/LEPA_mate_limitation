"""breeding_system_cartoon.py — what we know about the breeding system of
slickspot peppergrass (Lepidium papilliferum, LEPA), for a broad biology
audience: self-incompatibility (SI) under genetic control, its escape
(self-compatibility, SC), and the cost of each.

Two rows — SI on top, SC below; mechanism on the left, impact on the right:
  1. **Self-incompatible breeding system: plants reject their own pollen** —
     the S-locus (a male gene, SCR, and a female gene, SRK) gives each plant an
     S-type; pollen carrying the flower's own S-type is rejected (no seed),
     pollen of another S-type is accepted (seed).
  2. **Fragmentation → genetic drift → mate limitation → population decline** —
     small, isolated populations lose S-types, plants run out of compatible
     mates and set little seed (mate limitation).
  3. **Some plants broke the gene** — plants whose S-locus copies are broken
     accept their own pollen and set seed alone (self-compatible).
  4. **Selfing → inbreeding → reduced fitness → population decline** —
     selfed offspring grow more alike and weaker
     generation after generation (inbreeding depression).

Flower, self-pollen loop and silique glyphs are borrowed from the ISI
breeding-system cartoon of SRK_bioinformatics (Canu_amplicon,
draw_ISI_cartoon.py); the gene-copy icon follows its SRK functional-copy
figure. LEPA numbers from SRK_bioinformatics (Canu-amplicon genotypes:
263 plants, 32 S-types). Other icons and colours as restoration_cartoon.py.
The challenge (restoration_program_cartoon.panel_challenge, shared with
Figure 2) sits on top: the problem these panels explain.
Output: figures_seed_bulking/breeding_system_cartoon.{png,pdf}, plus talk slides:
_slide1 (the challenge alone, 16:9), _slide2 (panels 1, 2) and _slide3
(panels 3, 4), with _slide2_step1 and _slide3_step1 showing the first panel.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch

import restoration_cartoon as rc
from restoration_cartoon import (ALLELE, MATE_LIM, OUT, OVER, SPOT, arrow, bee, cross_mark, dna,
                                 frame, lepa)

# palette of the ISI breeding-system figure (SRK_bioinformatics)
COL_PETAL, COL_PETAL_EDGE = "#FFFFFF", "#B8B8B8"
COL_CENTER, COL_ANTHER_DOT = "#F4C430", "#6B4E1E"
COL_STEM, COL_LEAF = "#4CAF50", "#66BB6A"
COL_POD, COL_POD_EDGE = "#F5EFE0", "#7D5A2B"
COL_SEED, COL_SEED_EDGE = "#A0522D", "#5D3A1A"
COL_NO = "#B0B0B0"
SI_RED, SC_GREEN = "#b2182b", "#1b7837"
SELF_PURPLE = "#5D4E7A"
GENE_ON, GENE_OFF = "#2166AC", "#555555"


# ---- glyphs borrowed from draw_ISI_cartoon.py (data units, 0–100 axes) ------
def flower(ax, cx, cy, r=6.0, stem_bottom=None, s_type=None):
    """4-petal white Brassicaceae flower; optional stem and S-type badge."""
    if stem_bottom is not None:
        ax.plot([cx, cx], [stem_bottom, cy - r * 0.55], color=COL_STEM, lw=4,
                solid_capstyle="round", zorder=1)
        ax.add_patch(Ellipse((cx - r * 0.75, (cy + stem_bottom) / 2), r * 1.0, r * 0.32, angle=30,
                             fc=COL_LEAF, ec=COL_STEM, lw=1, zorder=2))
    for angle in (0, 90, 180, 270):
        a = np.deg2rad(angle)
        ax.add_patch(Ellipse((cx + r * 0.55 * np.cos(a), cy + r * 0.55 * np.sin(a)), r * 0.85,
                             r * 0.55, angle=angle, fc=COL_PETAL, ec=COL_PETAL_EDGE, lw=1.5,
                             zorder=3))
    ax.add_patch(Circle((cx, cy), r * 0.24, fc=COL_CENTER, ec="#B8930A", lw=1.2, zorder=4))
    for a in np.deg2rad((45, 135, 225, 315)):
        ax.add_patch(Circle((cx + r * 0.13 * np.cos(a), cy + r * 0.13 * np.sin(a)), r * 0.05,
                            fc=COL_ANTHER_DOT, ec="none", zorder=5))
    if s_type is not None:
        s_badge(ax, cx + r * 0.95, cy - r * 0.95, s_type)


def silique_full(ax, cx, cy, size=5.0):
    ax.add_patch(Ellipse((cx, cy), size * 1.3, size * 2.1, fc=COL_POD, ec=COL_POD_EDGE, lw=2,
                         zorder=2))
    ax.plot([cx, cx], [cy - size * 0.9, cy + size * 0.9], color=COL_POD_EDGE, lw=1,
            ls=(0, (2, 2)), zorder=2)
    for py in (0.65, 0.30, -0.05, -0.40, -0.72):
        for dx in (-0.22, 0.22):
            ax.add_patch(Ellipse((cx + dx * size, cy + py * size), size * 0.28, size * 0.20,
                                 fc=COL_SEED, ec=COL_SEED_EDGE, lw=1, zorder=3))


def silique_empty(ax, cx, cy, size=5.0):
    ax.add_patch(Ellipse((cx, cy), size * 1.3, size * 2.1, fc="none", ec=COL_NO, lw=2,
                         ls=(0, (5, 3)), zorder=2))
    ox, oy = size * 0.55, size * 0.9
    for sgn in (1, -1):
        ax.plot([cx - ox, cx + ox], [cy - sgn * oy, cy + sgn * oy], color=SI_RED, lw=4,
                solid_capstyle="round", zorder=3)


def self_loop(ax, cx, cy, r, label=True):
    y = cy + r * 0.4
    ax.add_patch(FancyArrowPatch((cx + r * 0.55, y), (cx - r * 0.55, y),
                                 connectionstyle="arc3,rad=1.3", arrowstyle="->",
                                 mutation_scale=22, color=SELF_PURPLE, lw=2.2, zorder=6))
    if label:
        ax.text(cx, cy + r * 1.45, "own pollen", ha="center", va="bottom", fontsize=12,
                color=SELF_PURPLE, style="italic", fontweight="bold")


# ---- new glyphs --------------------------------------------------------------
def s_badge(ax, x, y, k, r=2.0):
    """A plant's S-type: colored disc (one color = one S-type)."""
    ax.add_patch(Circle((x, y), r, fc=ALLELE[k], ec="k", lw=0.8, zorder=7))


def gene_copies(ax, x, y, working, w=4.2, h=11):
    """Four copies of the S-locus (LEPA has four copies of each chromosome):
    blue = working, grey with red × = broken (after SRK_bioinformatics)."""
    for k, ok in enumerate(working):
        x0 = x + k * (w + 1.2) + (1.6 if k > 1 else 0)
        ax.add_patch(FancyBboxPatch((x0, y), w, h, boxstyle="round,pad=0,rounding_size=1.2",
                                    fc="#D6E0EA" if k < 2 else "#E3DAE6", ec="#555", lw=1,
                                    zorder=3))
        cx, cy = x0 + w / 2, y + h * 0.62
        ax.add_patch(Ellipse((cx, cy), w * 0.8, h * 0.26, fc=GENE_ON if ok else GENE_OFF,
                             ec="k", lw=0.9, zorder=4))
        if not ok:
            for sgn in (1, -1):
                ax.plot([cx - w * 0.45, cx + w * 0.45], [cy - sgn * h * 0.16, cy + sgn * h * 0.16],
                        color=SI_RED, lw=2, zorder=5)


def tag(ax, x, y, text, col):
    ax.text(x, y, text, ha="center", va="center", fontsize=13, fontweight="bold", color="white",
            bbox=dict(boxstyle="round,pad=0.35", fc=col, ec="none"), zorder=10)


# ---- panels ------------------------------------------------------------------
def panel_si(ax):
    frame(ax, 1, "Self-incompatible breeding system:\n    plants reject their own pollen",
          "the S-locus (a male + a female gene) gives each plant an\n"
          "S-type: pollen with the flower's own S-type is rejected")
    tag(ax, 82, 93, "self-incompatible", SI_RED)
    dna(ax, 6, 87.5, h=10.5, w=2.0)                      # S-types are under genetic control
    ax.text(3.5, 84.5, "under genetic control", ha="left", va="center", fontsize=11, color="#333",
            style="italic")
    arrow(ax, (9.5, 93), (13, 93), col="#666", lw=1.6)
    for k, t in enumerate((0, 2, 4)):                    # key: colour = S-type
        s_badge(ax, 16 + 4.5 * k, 93, t, r=1.6)
    ax.text(29, 93, "color = S-type\n(32 known in LEPA)", ha="left", va="center", fontsize=11,
            linespacing=1.1, color="#333")
    # own pollen → no seed
    flower(ax, 24, 66, r=7, stem_bottom=50, s_type=0)
    self_loop(ax, 24, 66, 7)
    arrow(ax, (24, 46), (24, 37), col="#666", lw=1.8)
    silique_empty(ax, 24, 26, size=5)
    ax.text(24, 15.5, "no seed", ha="center", va="top", fontsize=13, color=SI_RED,
            fontweight="bold")
    # pollen of another S-type → seed
    flower(ax, 58, 66, r=6, stem_bottom=52, s_type=2)
    flower(ax, 84, 66, r=7, stem_bottom=50, s_type=0)
    ax.add_patch(FancyArrowPatch((61, 70), (81, 70), connectionstyle="arc3,rad=-0.5",
                                 arrowstyle="->", mutation_scale=20, color="#2C5F8D", lw=2,
                                 zorder=6))
    ax.text(71, 80, "other S-type", ha="center", va="bottom", fontsize=12, color="#2C5F8D",
            style="italic", fontweight="bold")
    arrow(ax, (84, 46), (84, 37), col="#666", lw=1.8)
    silique_full(ax, 84, 26, size=5)
    ax.text(84, 15.5, "seed", ha="center", va="top", fontsize=13, color=COL_POD_EDGE,
            fontweight="bold")
    ax.plot([42, 42], [18, 86], color="#DDD", lw=1.2, zorder=0)


def panel_mates(ax):
    frame(ax, 2, "Fragmentation → genetic drift\n    → mate limitation → population decline",
          "small, isolated populations lose S-types by drift:\n"
          "plants run out of compatible mates and set little seed")
    tag(ax, 79, 93, "mate limitation", MATE_LIM)
    for x0, ok, types in ((26, True, (0, 2, 4, 6, 1)), (74, False, (0, 0, 1, 0, 1))):
        ax.add_patch(Ellipse((x0, 58), 42, 30, fc=SPOT, ec="#B08D57", lw=1.5, zorder=1))
        ax.text(x0, 77, "many S-types" if ok else "few S-types", ha="center", fontsize=13,
                fontweight="bold", color="#333")
        xs = np.linspace(x0 - 14, x0 + 14, len(types))
        for k, (xx, t) in enumerate(zip(xs, types)):
            yy = 52 + (5 if k % 2 else 0)
            lepa(ax, xx, yy, h=9, fruits=ok)
            s_badge(ax, xx, yy - 2.2, t, r=1.6)
        for a, b in ((0, 1), (1, 2), (2, 3), (3, 4)):          # pollinator visits
            (x1, x2), y1 = (xs[a], xs[b]), 67
            ax.add_patch(FancyArrowPatch((x1, y1), (x2, y1), arrowstyle="-", color="#666",
                                         lw=1.1, ls=(0, (3, 2)), connectionstyle="arc3,rad=-0.5",
                                         zorder=6))
        bee(ax, x0 + 2, 72, s=1.6)
        for k in range(3):
            (silique_full if ok else silique_empty)(ax, x0 - 9 + 9 * k, 28, size=3.4)
    arrow(ax, (41, 82), (59, 82), col="#555", lw=2.2)    # many → few S-types
    ax.text(50, 85, "habitat fragmentation\n→ genetic drift", ha="center", va="bottom",
            fontsize=11, color="#333", fontweight="bold", linespacing=1.05)
    cross_mark(ax, 74 - 7, 70.5, s=1.4)
    cross_mark(ax, 74 + 7, 70.5, s=1.4)
    ax.text(74, 17, "in the field: less seed than\nplant size predicts", ha="center", va="center",
            fontsize=11.5, color=MATE_LIM, style="italic", linespacing=1.1)


def panel_sc(ax):
    frame(ax, 3, "Rapid adaptation: escape from SI\n    → self-compatible, short-term gain",
          "plants with broken S-locus copies accept their own pollen\n"
          "and set seed alone, without mates")
    tag(ax, 82, 93, "self-compatible", SC_GREEN)
    gene_copies(ax, 6, 64, [True] * 4)
    ax.text(17.5, 60, "working S-locus\n(4 copies)", ha="center", va="top", fontsize=11.5,
            linespacing=1.1)
    arrow(ax, (31, 69.5), (37, 69.5), col="#666", lw=1.8)
    gene_copies(ax, 39, 64, [False] * 4)
    ax.text(50.5, 60, "broken copies", ha="center", va="top", fontsize=11.5, color=SI_RED)
    flower(ax, 80, 66, r=7, stem_bottom=50)
    self_loop(ax, 80, 66, 7)
    arrow(ax, (80, 46), (80, 37), col="#666", lw=1.8)
    silique_full(ax, 76, 26, size=5)                     # more seed: short-term gain
    silique_full(ax, 84, 26, size=5)
    ax.text(80, 15.5, "seed ↑", ha="center", va="top", fontsize=14.5, color=SC_GREEN,
            fontweight="bold")
    # LEPA numbers (SRK_bioinformatics, 263 genotyped plants)
    ax.text(30, 44, "LEPA — 263 plants genotyped", ha="center", fontsize=12, fontweight="bold",
            color="#333")
    for k, (n, lab, col) in enumerate(((247, "self-incompatible", SI_RED),
                                       (15, "partly broken", "#E08214"),
                                       (1, "self-compatible", SC_GREEN))):
        y = 37 - 7 * k
        ax.text(17, y, f"{n}", ha="right", va="center", fontsize=13, fontweight="bold", color=col)
        ax.text(19, y, lab, ha="left", va="center", fontsize=12, color=col)
    ax.text(30, 16, "in the field: more seed than\nplant size predicts", ha="center", va="center",
            fontsize=11.5, color=OVER, style="italic", linespacing=1.1)


def panel_inbreeding(ax):
    frame(ax, 4, "Long-term cost: inbreeding depression\n    → population decline",
          "selfed offspring grow more alike and weaker,\n"
          "generation after generation")
    tag(ax, 82, 93, "inbreeding", OVER)
    gens = [(17, 26, (0, 2), True, "parent"), (50, 20, (0, 0), True, "selfed\noffspring"),
            (82, 14, (0, 0), False, "next\ngeneration")]
    for x, h, types, fr, lab in gens:
        lepa(ax, x, 42, h=h, fruits=fr, faded=not fr)
        for k, t in enumerate(types):
            s_badge(ax, x + (k - 0.5) * 4.6, 38, t, r=2.0)
        ax.text(x, 33, lab, ha="center", va="top", fontsize=12, linespacing=1.1)
    for xa, xb in ((29, 40), (61, 72)):
        arrow(ax, (xa, 70), (xb, 66), col=SELF_PURPLE, lw=2.4, rad=-0.4)
        ax.text((xa + xb) / 2, 76, "self", ha="center", fontsize=12, color=SELF_PURPLE,
                style="italic", fontweight="bold")
    for k, lab in enumerate(("seed ↓", "survival ↓", "diversity ↓")):
        ax.text(17 + 32.5 * k, 18, lab, ha="center", va="center", fontsize=13.5,
                fontweight="bold", color=OVER)


PANELS = (panel_si, panel_mates, panel_sc, panel_inbreeding)


def draw():
    from restoration_program_cartoon import panel_challenge      # (lazy: avoids an import cycle)
    rc.RNG = np.random.default_rng(5)
    fig = plt.figure(figsize=(16, 19.0))
    top = fig.add_axes([0.02, 0.825, 0.96, 0.96 * 16 / 19.0 * 38 / 200])
    panel_challenge(top, "The challenge: very little seed to start with — why?")
    gs = fig.add_gridspec(2, 2, left=0.02, right=0.98, top=0.735, bottom=0.015, wspace=0.1,
                          hspace=0.17)
    axes = np.array([[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(2)])
    for ax, f in zip(axes.flat, PANELS):
        f(ax)
    fig.canvas.draw()
    to_fig = lambda ax, xy: fig.transFigure.inverted().transform(ax.transData.transform(xy))
    style = dict(transform=fig.transFigure, arrowstyle="-|>", mutation_scale=26, color="#333",
                 lw=3)
    p, q = to_fig(top, (100, 0)), to_fig(axes[0, 0], (50, 100))  # challenge → why
    fig.add_artist(FancyArrowPatch((p[0], p[1] - 0.004), (p[0], q[1] + 0.03), **style))
    for a, b in axes:                                    # mechanism → impact, per row
        fig.add_artist(FancyArrowPatch(to_fig(a, (100, 50)), to_fig(b, (0, 50)), **style))
    return fig


def draw_challenge_slide():
    """16:9 slide: the challenge alone, edge to edge."""
    from restoration_program_cartoon import panel_challenge_slide
    rc.RNG = np.random.default_rng(5)
    fig = plt.figure(figsize=(16, 9))
    panel_challenge_slide(fig.add_axes([0, 0, 1, 1]))
    return fig


def draw_pair_slide(left, right=None):
    """Landscape slide: two square panels, cropped close to the titles and
    panels (a slim white edge, for a colored slide background); the right
    panel blank (same layout) until it is shown."""
    rc.RNG = np.random.default_rng(5)
    side, gap, title, m = 7.6, 0.5, 0.68, 0.12           # inches
    x0, h = (16 - 2 * side - gap) / 2, m + side + title + m
    fig = plt.figure(figsize=(16, h))
    a, b = (fig.add_axes([x / 16, m / h, side / 16, side / h]) for x in (x0, x0 + side + gap))
    left(a)
    if right is None:
        b.axis("off")
        return fig
    right(b)
    fig.canvas.draw()
    to_fig = lambda ax, xy: fig.transFigure.inverted().transform(ax.transData.transform(xy))
    fig.add_artist(FancyArrowPatch(to_fig(a, (100, 50)), to_fig(b, (0, 50)),
                                   transform=fig.transFigure, arrowstyle="-|>", mutation_scale=26,
                                   color="#333", lw=3))
    return fig


def main():
    OUT.mkdir(exist_ok=True)
    fig = draw()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"breeding_system_cartoon.{ext}", dpi=200 if ext == "png" else None)
    plt.close(fig)
    slides = {"slide1": draw_challenge_slide,
              "slide2_step1": lambda: draw_pair_slide(panel_si),
              "slide2": lambda: draw_pair_slide(panel_si, panel_mates),
              "slide3_step1": lambda: draw_pair_slide(panel_sc),
              "slide3": lambda: draw_pair_slide(panel_sc, panel_inbreeding)}
    for name, make in slides.items():
        fig = make()
        fig.savefig(OUT / f"breeding_system_cartoon_{name}.png", dpi=200)
        plt.close(fig)
    print(f"  wrote {OUT}/breeding_system_cartoon.png (+ talk slides _slide1–3)")


if __name__ == "__main__":
    main()
