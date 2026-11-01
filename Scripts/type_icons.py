"""type_icons.py — small badge icons for the three reproductive types, to
help read the seed-set figures (seed by size, population seed yield).

  Mate-limited (magenta)      two flowers of the same S-type: pollen
                              rejected (✗), empty fruit — no compatible mate.
  Reproductively intact (teal) two flowers of different S-types: pollen
                              accepted (✓), full fruit — seed as expected.
  Over-producing (gold)       one flower pollinating itself (own-pollen
                              loop), extra-full fruit — selfing, more seed.

Drawn with the flower, pollen-arc and fruit glyphs of
breeding_system_cartoon.py (after the ISI figures of SRK_bioinformatics),
inside a round badge ringed in the type's colour. Use add_icon(fig, kind,
[x, y, w, h]) to place one on any figure (figure coordinates).
Output (preview sheet): figures_seed_bulking/type_icons.png
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

from breeding_system_cartoon import (ALLELE, COL_POD_EDGE, SELF_PURPLE, flower, silique_empty,
                                     silique_full)
from restoration_cartoon import INTACT, MATE_LIM, OUT, OVER

KINDS = {"Mate-limited": MATE_LIM, "Reproductively intact": INTACT, "Over-producing": OVER}


def _scale(ax):
    """Line widths and font sizes follow the icon's size (1 = a 2-inch icon)."""
    return ax.get_position().height * ax.figure.get_figheight() / 2.0


def _badge(ax, col, k):
    ax.set_xlim(0, 100), ax.set_ylim(0, 100), ax.set_aspect("equal"), ax.axis("off")
    ax.add_patch(Circle((50, 50), 47, fc="white", ec=col, lw=4 * k, zorder=0))


def _two_flowers(ax, s1, s2, ok, col, k):
    for x, s in ((29, s1), (71, s2)):
        flower(ax, x, 60, r=19)
        ax.add_patch(Circle((x + 13, 47), 6.5, fc=ALLELE[s], ec="k", lw=0.8 * k, zorder=7))
    ax.add_patch(FancyArrowPatch((35, 73), (65, 73), connectionstyle="arc3,rad=-0.55",
                                 arrowstyle="->", mutation_scale=16 * k, color=col, lw=2.6 * k,
                                 ls="-" if ok else (0, (2.5, 1.5)), zorder=6))
    ax.text(50, 83, "✓" if ok else "✗", ha="center", va="center", fontsize=22 * k,
            fontweight="bold", color=col, zorder=8,
            bbox=dict(boxstyle="circle,pad=0.05", fc="white", ec="none"))


def icon(ax, kind):
    col, k = KINDS[kind], _scale(ax)
    _badge(ax, col, k)
    if kind == "Mate-limited":
        _two_flowers(ax, 0, 0, False, col, k)             # same S-type: rejected
        silique_empty(ax, 50, 22, size=11)
    elif kind == "Reproductively intact":
        _two_flowers(ax, 0, 2, True, col, k)              # different S-types: accepted
        silique_full(ax, 50, 22, size=11)
    else:                                                  # over-producing: selfing
        flower(ax, 36, 58, r=22)
        ax.add_patch(FancyArrowPatch((46, 70), (26, 70), connectionstyle="arc3,rad=1.3",
                                     arrowstyle="->", mutation_scale=16 * k, color=SELF_PURPLE,
                                     lw=2.6 * k, zorder=6))
        for dx in (-8, 8):
            silique_full(ax, 72 + dx, 42, size=10.5)
        ax.text(72, 18, "×2", ha="center", va="center", fontsize=17 * k, fontweight="bold",
                color=COL_POD_EDGE)


def add_icon(fig, kind, rect):
    """Badge icon for `kind` in figure coordinates rect = [x, y, w, h]."""
    ax = fig.add_axes(rect)
    icon(ax, kind)
    return ax


def main():
    fig = plt.figure(figsize=(9, 3.6))
    for k, kind in enumerate(KINDS):
        add_icon(fig, kind, [0.03 + 0.33 * k, 0.2, 0.27, 0.75])
        fig.text(0.165 + 0.33 * k, 0.1, kind.replace("Reproductively ", "Reproductively\n"),
                 ha="center", va="center", fontsize=12, fontweight="bold", color=KINDS[kind])
    fig.savefig(OUT / "type_icons.png", dpi=220)
    plt.close(fig)
    print(f"  wrote {OUT}/type_icons.png")


if __name__ == "__main__":
    main()
