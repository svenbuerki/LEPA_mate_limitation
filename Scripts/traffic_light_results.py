"""traffic_light_results.py — the traffic light of Figure 1 filled with our
results: every population with two seasons of seed data (2025, 2026) under
its light, as a chip with its EO —

  red     AVOID    over-producing populations
  orange  TARGET   mate-limited populations (hollow chip = candidate: below
                   expectation, but not significantly in either season)
  green   INTACT   reproductively intact populations (the seed sources)

and, in the gray footer, the populations not yet under a light: switching
between seasons (not yet clear) and one season of data so far, with the
candidate sources whose 2026 seed is being cleaned (one season, 2025 seed set
significantly within half – twice, 2026 mothers collected — the rule of the
status map, step32_seed_bulking_maps.R). Read from
tables/step32_seed_bulking_populations.tsv (step32_seed_bulking.py), so it
follows the data. Colors, badges and action symbols as Figures 1 – 4.
Output: figures_seed_bulking/traffic_light_results.{png,pdf} (16:9, also the
talk slide). Run from the repository root:
    python Scripts/traffic_light_results.py
"""
from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Circle, FancyBboxPatch, Polygon

from restoration_cartoon import INTACT, MATE_LIM, OUT, OVER
from restoration_program_cartoon import shield
from type_icons import add_icon, lifebuoy, no_sign

TABLE = "tables/step32_seed_bulking_populations.tsv"
BL_TABLE = "tables/step30h_bl_definition.tsv"           # populations per BL: all 44
W, H = 16.0, 9.0                                         # inches; data units 0–160 × 0–90
ROWS = (  # lamp color, type, category, biological name, y
    (OVER, "Over-producing", "AVOID", "over-producing", 67),
    (MATE_LIM, "Mate-limited", "TARGET", "mate-limited", 45),
    (INTACT, "Reproductively intact", "INTACT", "reproductively intact", 23),
)
CHIP_W, CHIP_H, CHIP_GAP, CHIP_X0 = 11.2, 10.5, 1.2, 74.0


def eo_label(codes: str) -> str:
    """EO of a population, without the part suffix (EO27-1 → EO27)."""
    eos = dict.fromkeys(re.sub(r"(EO\d+).*", r"\1", c) for c in str(codes).split(","))
    return " + ".join(eos)


def chip(ax, x, y, pid, eo, bl, col, hollow):
    ax.add_patch(FancyBboxPatch((x, y - CHIP_H / 2), CHIP_W, CHIP_H,
                                boxstyle="round,pad=0,rounding_size=1.6",
                                fc="white" if hollow else col, ec=col, lw=2.2, zorder=3))
    fg = col if hollow else "white"
    ax.text(x + CHIP_W / 2, y + 2.3, f"P{pid}", ha="center", va="center", fontsize=15.5,
            fontweight="bold", color=fg, zorder=4)
    ax.text(x + CHIP_W / 2, y - 2.7, eo, ha="center", va="center", fontsize=14,
            color=fg, zorder=4)


def draw() -> None:
    st = pd.read_csv(TABLE, sep="\t")
    st["candidate"] = st["candidate"].astype(str).isin(["True", "true", "1"])
    one = st[st["type"] == "One season only"]
    cand_src = one[(one["n_plants_2026"] == 0) & (one["qeq_2025"] < 0.05)
                   & (one["mothers_collected_2026"] > 0)]
    n_two = int((st["type"] != "One season only").sum())
    n_all = int(pd.read_csv(BL_TABLE, sep="\t")["n_populations"].sum())
    n_one = n_all - n_two                                # incl. those seen without seed data yet

    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 160), ax.set_ylim(0, 90), ax.set_aspect("equal"), ax.axis("off")
    ax.text(80, 86.5, "Our populations under the traffic light", ha="center", va="center",
            fontsize=27, fontweight="bold")
    ax.text(80, 81.3, f"the {n_two} populations with two seasons of seed set for plant size "
            "(2025, 2026)", ha="center", va="center", fontsize=15.5, color="#444",
            style="italic")
    # the traffic light of Figure 1
    ax.add_patch(FancyBboxPatch((16.3, 1), 3.4, 10, boxstyle="square,pad=0", fc="#455A64",
                                ec="none", zorder=1))
    ax.add_patch(FancyBboxPatch((7, 12), 22, 66, boxstyle="round,pad=0,rounding_size=3.5",
                                fc="#263238", ec="#111", lw=2, zorder=2))
    to_fig = (ax.transData + fig.transFigure.inverted()).transform
    for col, kind, cat, name, y in ROWS:
        ax.add_patch(Polygon([[9, y + 8.6], [27, y + 8.6], [28.5, y + 10], [7.5, y + 10]],
                             fc="#111", ec="none", zorder=3))
        ax.add_patch(Circle((18, y), 8.0, fc=col, ec="#111", lw=2, zorder=4))
        ax.add_patch(Circle((15.4, y + 2.7), 2.1, fc="white", alpha=0.35, ec="none", zorder=5))
        if y != 23:
            ax.plot([33, 154], [y - 11] * 2, color="#DDD", lw=1, zorder=0)
        (x0, y0), (x1, y1) = to_fig([(33, y - 6.5), (46, y + 6.5)])
        add_icon(fig, kind, [x0, y0, x1 - x0, y1 - y0])
        pops = st[st["type"] == kind].sort_values(["candidate", "populationID"])
        ax.text(49, y + 5, cat, ha="left", va="center", fontsize=24, fontweight="bold",
                color=col)
        ax.text(49, y - 0.5, name, ha="left", va="center", fontsize=12, color=col,
                fontweight="bold")
        n = len(pops)
        ax.text(49, y - 5.5, f"{n} population{'s' if n != 1 else ''}", ha="left",
                va="center", fontsize=13, color="#333")
        for k, (_, r) in enumerate(pops.iterrows()):
            chip(ax, CHIP_X0 + k * (CHIP_W + CHIP_GAP), y, int(r["populationID"]),
                 eo_label(r["locationCodes"]), r["BL"], col, bool(r["candidate"]))
    no_sign(ax, 153.5, 67, r=4.8)                        # same symbols as Figures 1 and 3
    lifebuoy(ax, 153.5, 45, r=4.8)
    shield(ax, 153.5, 23.5, s=0.54)
    # populations not yet under a light
    unclear = st[st["type"] == "Not yet clear"]
    src = ", ".join(f"P{p}" for p in sorted(cand_src["populationID"]))
    lines = [f"Not yet under a light: {n_one} populations with one season of data so far "
             f"(likely sources being confirmed: {src})"]
    if len(unclear):
        lines.append("and " + ", ".join(f"P{int(p)} ({eo_label(c)})" for p, c in
                                        zip(unclear["populationID"], unclear["locationCodes"]))
                     + " — switches between seasons · hollow chip = candidate "
                     "(below expectation, not yet significant)")
    ax.text(93, 4.8, "\n".join(lines), ha="center", va="center", fontsize=11.5, color="#333",
            linespacing=1.35,
            bbox=dict(boxstyle="round,pad=0.5", fc="#F3F4F6", ec="#999", lw=1))
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"traffic_light_results.{ext}", dpi=200)
    plt.close(fig)
    print(f"  wrote {OUT}/traffic_light_results.{{png,pdf}}")


if __name__ == "__main__":
    draw()
