"""traffic_light_cartoon.py — Figure 1 of Seed_bulking.md: the traffic light used
to prioritize EOs for the species recovery plan (Recovery Criterion 2; USFWS
2026). A literal traffic light, one lamp per category, in plain language for
every audience (land managers first); the biology behind each category
(mate limitation, inbreeding depression) comes later, with Figure 2:

  red     AVOID    do not use as a seed source (plants set more seed than
                   expected, likely by pollinating themselves)
  orange  TARGET   restore: bring in new genetic diversity (plants set less
                   seed than expected)
  green   INTACT   conserve: the seed source for restoration (plants set the
                   seed expected for their size)

Colors and action symbols (prohibition sign, ring buoy, shield) as on every
other figure. Output: figures_seed_bulking/traffic_light_cartoon.{png,pdf}
(16:9, also the talk slide). Run from the repository root:
    python Scripts/traffic_light_cartoon.py
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, Polygon

from restoration_cartoon import INTACT, MATE_LIM, OUT, OVER
from restoration_program_cartoon import shield
from type_icons import lifebuoy, no_sign

W, H = 16.0, 9.0                                         # inches; data units 0–160 × 0–90
ROWS = (  # lamp color, category, action, what we see, y
    (OVER, "AVOID", "do not use as a seed source",
     "plants set more seed than expected for their size,\nlikely by pollinating themselves", 68),
    (MATE_LIM, "TARGET", "restore — bring in new genetic diversity",
     "plants set less seed than expected for their size", 45),
    (INTACT, "INTACT", "conserve — the seed source for restoration",
     "plants set the seed expected for their size", 22),
)


def draw() -> None:
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 160), ax.set_ylim(0, 90), ax.set_aspect("equal"), ax.axis("off")
    ax.text(80, 86.5, "Prioritizing EOs for the species recovery plan",
            ha="center", va="center", fontsize=27, fontweight="bold")
    ax.text(80, 81.3, "Recovery Criterion 2: at least 100 EOs persist over 12 years, "
            "at least 75 of them in fair to excellent condition",
            ha="center", va="center", fontsize=15.5, color="#444", style="italic")
    # the traffic light: pole, housing, three lamps under their visors
    ax.add_patch(FancyBboxPatch((16.3, 1), 3.4, 10, boxstyle="square,pad=0", fc="#455A64",
                                ec="none", zorder=1))
    ax.add_patch(FancyBboxPatch((7, 11), 22, 67, boxstyle="round,pad=0,rounding_size=3.5",
                                fc="#263238", ec="#111", lw=2, zorder=2))
    for col, _, _, _, y in ROWS:
        ax.add_patch(Polygon([[9, y + 8.6], [27, y + 8.6], [28.5, y + 10], [7.5, y + 10]],
                             fc="#111", ec="none", zorder=3))                  # visor
        ax.add_patch(Circle((18, y), 8.0, fc=col, ec="#111", lw=2, zorder=4))
        ax.add_patch(Circle((15.4, y + 2.7), 2.1, fc="white", alpha=0.35, ec="none", zorder=5))
    # one row per lamp: category, action, what we see in the field, action symbol
    for col, cat, action, seen, y in ROWS:
        if y != 22:
            ax.plot([35, 154], [y - 11.5] * 2, color="#DDD", lw=1, zorder=0)
        ax.text(38, y + 5.8, cat, ha="left", va="center", fontsize=32, fontweight="bold",
                color=col)
        ax.text(38, y - 0.6, action, ha="left", va="center", fontsize=21, color=col,
                fontweight="bold")
        ax.text(38, y - 6.8, seen, ha="left", va="center", fontsize=15, color="#333",
                linespacing=1.2)
    shield(ax, 146, 23.5, s=0.62)                        # same symbols as Figures 3 and 4
    lifebuoy(ax, 146, 45, r=5.6)
    no_sign(ax, 146, 68, r=5.6)
    ax.text(93, 4.5, "Each EO is categorized by how its populations reproduce "
            "(seed set for plant size), then confirmed with genetics",
            ha="center", va="center", fontsize=12.5, color="#333",
            bbox=dict(boxstyle="round,pad=0.5", fc="#F3F4F6", ec="#999", lw=1))
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"traffic_light_cartoon.{ext}", dpi=200)
    plt.close(fig)
    print(f"  wrote {OUT}/traffic_light_cartoon.{{png,pdf}}")


if __name__ == "__main__":
    draw()
