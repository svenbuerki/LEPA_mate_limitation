"""Step 09 — Bottleneck Lineage maps on the Snake River Plain.

1. Overview: one dot per population, coloured by BL, sized by fertile
   plants (2025 + 2026), population IDs labelled, BL convex hulls.
2. One map per BL: EO convex hulls and population polygons on the BL
   extent, plus one zoom panel per population showing its 50 m demes in
   each season (2025 dashed / open markers, 2026 solid / filled).
   Polygons are convex hulls of the events, buffered in metres (EOs
   50 m, populations and demes 25 m = half the pollinator radius) so
   single-event units remain visible.

Maps carry coordinates (figures only; no coordinate tables are
published).

Outputs
-------
figures/step30h_overview_by_BL.png / .pdf
figures/step30h_BL_map_BL{1..5}.png, figures/step30h_BL_maps.pdf (all BLs)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial import ConvexHull

from srk_bl_constants import BL_COLORS, BL_ORDER, base_eo
from srk_units import load_event_units

TABLES = Path("tables")
FIGURES = Path("figures")
CROSSWALK_TSV = TABLES / "step29a_population_crosswalk.tsv"
EO_BUFFER_M = 50.0
UNIT_BUFFER_M = 25.0
M_PER_DEG_LAT = 110_540.0


def main() -> None:
    pts = pd.read_csv(TABLES / "step30h_cluster_assignments.tsv",
                       sep="\t", encoding="utf-8-sig")
    bl_def = pd.read_csv(TABLES / "step30h_bl_definition.tsv",
                           sep="\t", encoding="utf-8-sig")
    # Phase III remapped cluster_assignments: it already carries the
    # new populationID + BL columns.

    fig, ax = plt.subplots(figsize=(11.5, 7.5))
    size_min, size_max = 30, 500
    n_max = max(pts["n_fert_total"].max(), 1)
    sizes = size_min + (size_max - size_min) * np.sqrt(
        pts["n_fert_total"] / n_max
    )

    # BL hulls
    for bl, grp in pts.groupby("BL"):
        color = BL_COLORS.get(bl, "#777777")
        if len(grp) >= 3:
            hull_pts = grp[["lon", "lat"]].to_numpy()
            try:
                hull = ConvexHull(hull_pts)
                poly = hull_pts[hull.vertices]
                poly = np.vstack([poly, poly[:1]])
                ax.fill(poly[:, 0], poly[:, 1], color=color,
                        alpha=0.08, zorder=1)
                ax.plot(poly[:, 0], poly[:, 1], color=color,
                        linewidth=1.3, linestyle="--", alpha=0.75,
                        zorder=2)
            except Exception:
                pass

    # Population dots coloured by BL
    for bl, grp in pts.groupby("BL"):
        color = BL_COLORS.get(bl, "#777777")
        s = sizes.loc[grp.index].to_numpy()
        row = bl_def[bl_def["BL"] == bl].iloc[0]
        ax.scatter(grp["lon"], grp["lat"], s=s, c=color,
                   alpha=0.85, edgecolor="white", linewidth=0.7,
                   zorder=3,
                   label=f"{bl} (n = {int(row['n_populations'])}, "
                         f"{row['convex_hull_area_km2']:.0f} km²)")

    # Label populations (only when there's room — skip the densest BL5)
    for _, r in pts.iterrows():
        ax.annotate(f"P{int(r['populationID'])}",
                    (r["lon"], r["lat"]),
                    xytext=(5, 3), textcoords="offset points",
                    fontsize=7, color="#333", alpha=0.85)

    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title(
        "Populations across the Snake River Plain, grouped into five "
        "Bottleneck Lineages (BL; Ward's clustering, k = 5)\n"
        "Dot size ∝ fertile plants counted in 2025 + 2026; "
        "dashed outlines = BL convex hulls.",
        fontsize=11,
    )
    ax.legend(loc="upper right", fontsize=9, frameon=True,
              title="Bottleneck Lineage")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_png = FIGURES / "step30h_overview_by_BL.png"
    out_pdf = FIGURES / "step30h_overview_by_BL.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30h-V] Wrote {out_png.name} + .pdf")
    plot_bl_maps()


# ---------------------------------------------------------------------------
# Per-BL maps: EOs, populations, demes
# ---------------------------------------------------------------------------
def _to_m(lon, lat, lat0):
    k = M_PER_DEG_LAT * np.cos(np.radians(lat0))
    return np.asarray(lon) * k, np.asarray(lat) * M_PER_DEG_LAT


def _to_deg(x, y, lat0):
    k = M_PER_DEG_LAT * np.cos(np.radians(lat0))
    return np.asarray(x) / k, np.asarray(y) / M_PER_DEG_LAT


def buffered_hull(lon, lat, buffer_m, lat0, n_circle=24):
    """Convex hull of the points, buffered by `buffer_m` metres
    (Minkowski sum with a circle); returns (lon, lat) of a closed ring."""
    x, y = _to_m(lon, lat, lat0)
    a = np.linspace(0, 2 * np.pi, n_circle, endpoint=False)
    px = (x[:, None] + buffer_m * np.cos(a)).ravel()
    py = (y[:, None] + buffer_m * np.sin(a)).ravel()
    pts = np.column_stack([px, py])
    hull = ConvexHull(pts)
    ring = pts[hull.vertices]
    ring = np.vstack([ring, ring[:1]])
    return _to_deg(ring[:, 0], ring[:, 1], lat0)


def _scale_bar(ax, lat0, frac=0.25):
    """Round-number scale bar in the lower-left corner."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    width_m = (x1 - x0) * M_PER_DEG_LAT * np.cos(np.radians(lat0)) * frac
    nice = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]
    L = max([v for v in nice if v <= width_m] or [nice[0]])
    dx = L / (M_PER_DEG_LAT * np.cos(np.radians(lat0)))
    xs = x0 + 0.05 * (x1 - x0)
    ys = y0 + 0.06 * (y1 - y0)
    ax.plot([xs, xs + dx], [ys, ys], color="#222", lw=2.5,
            solid_capstyle="butt", zorder=6)
    ax.text(xs + dx / 2, ys + 0.02 * (y1 - y0),
            f"{L / 1000:g} km" if L >= 1000 else f"{L} m",
            ha="center", va="bottom", fontsize=7, zorder=6)


def _set_equal(ax, lat0):
    ax.set_aspect(1 / np.cos(np.radians(lat0)))


def plot_bl_maps() -> None:
    cw = pd.read_csv(CROSSWALK_TSV, sep="\t", encoding="utf-8-sig")
    units = load_event_units()
    ev = units.merge(cw[["event_year", "eventID", "lat", "lon"]],
                     on=["event_year", "eventID"], how="left")
    ev = ev.dropna(subset=["lat", "lon"])
    ev["EO"] = ev["locationCode"].map(
        lambda c: base_eo(c).replace("EO0", "EO") if base_eo(c) else c)
    lat0 = float(ev["lat"].mean())
    s_max = max(ev["n_fertile"].max(), 1)
    pdf = PdfPages(FIGURES / "step30h_BL_maps.pdf")
    for bl in BL_ORDER:
        sub = ev[ev["BL"] == bl]
        if sub.empty:
            continue
        col = BL_COLORS.get(bl, "#555")
        pops = sorted(sub["populationID"].unique())
        n_cols = 4
        n_rows = int(np.ceil(len(pops) / n_cols))
        # Overview height follows the BL's real shape (equal-area axes)
        bx, by = _to_m(sub["lon"], sub["lat"], lat0)
        aspect = (np.ptp(by) + 800) / (np.ptp(bx) + 800)
        over_h = float(np.clip(14 * aspect, 3.5, 9.0))
        fig = plt.figure(figsize=(16, over_h + 3.6 * n_rows + 0.8))
        gs = fig.add_gridspec(1 + n_rows, n_cols,
                              height_ratios=[over_h / 3.6] + [1] * n_rows,
                              hspace=0.35, wspace=0.25)
        ax = fig.add_subplot(gs[0, :])
        # EO hulls: every event of each EO present in this BL
        for eo in sorted(sub["EO"].unique()):
            e = ev[ev["EO"] == eo]
            lo, la = buffered_hull(e["lon"].to_numpy(), e["lat"].to_numpy(),
                                   EO_BUFFER_M, lat0)
            ax.fill(lo, la, color="#9a9a9a", alpha=0.10, zorder=1)
            ax.plot(lo, la, color="#6e6e6e", lw=1.0, ls="--", zorder=2)
            ax.annotate(eo, (lo.max(), la.max()), fontsize=8, color="#444",
                        xytext=(2, 2), textcoords="offset points", zorder=5)
        for p in pops:
            e = sub[sub["populationID"] == p]
            lo, la = buffered_hull(e["lon"].to_numpy(), e["lat"].to_numpy(),
                                   UNIT_BUFFER_M, lat0)
            ax.fill(lo, la, color=col, alpha=0.35, zorder=3)
            ax.plot(lo, la, color=col, lw=1.2, zorder=3)
            ax.annotate(f"P{p}", (e["lon"].mean(), e["lat"].mean()),
                        fontsize=8, weight="bold", color="#111",
                        xytext=(4, -10), textcoords="offset points", zorder=5)
        ax.margins(0.06)
        _set_equal(ax, lat0)
        _scale_bar(ax, lat0)
        ax.set_xlabel("Longitude (°)")
        ax.set_ylabel("Latitude (°)")
        ax.set_title(f"{bl} — element occurrences (grey dashed convex hulls) "
                     f"and populations (filled, P{{N}}); {len(pops)} populations",
                     loc="left", fontsize=11, weight="bold")
        ax.grid(True, alpha=0.25)

        # One zoom panel per population: demes per season
        for i, p in enumerate(pops):
            axp = fig.add_subplot(gs[1 + i // n_cols, i % n_cols])
            e = sub[sub["populationID"] == p]
            lo, la = buffered_hull(e["lon"].to_numpy(), e["lat"].to_numpy(),
                                   UNIT_BUFFER_M, lat0)
            axp.fill(lo, la, color=col, alpha=0.10, zorder=1)
            # deme outlines stay visible in large populations: buffer at
            # least 3 % of the panel extent
            ex, ey_ = _to_m(lo, la, lat0)
            deme_buf = max(UNIT_BUFFER_M, 0.03 * max(np.ptp(ex), np.ptp(ey_)))
            n_demes = {}
            for yr, ls, mfc in ((2025, "--", "white"), (2026, "-", col)):
                ey = e[e["event_year"] == yr]
                n_demes[yr] = ey["demeID"].nunique()
                for _, d in ey.groupby("demeID"):
                    dlo, dla = buffered_hull(d["lon"].to_numpy(),
                                             d["lat"].to_numpy(),
                                             deme_buf, lat0)
                    axp.plot(dlo, dla, color=col, lw=1.1, ls=ls, zorder=2)
                axp.scatter(ey["lon"], ey["lat"],
                            s=8 + 60 * np.sqrt(ey["n_fertile"] / s_max),
                            facecolors=mfc, edgecolors=col, linewidths=0.8,
                            zorder=3)
            axp.margins(0.15)
            _set_equal(axp, lat0)
            _scale_bar(axp, lat0, frac=0.35)
            codes = ", ".join(sorted(e["locationCode"].unique()))
            axp.set_title(f"P{p} {codes}\ndemes: 2025 {n_demes[2025]} · "
                          f"2026 {n_demes[2026]}", fontsize=8)
            axp.set_xticks([])
            axp.set_yticks([])
        fig.suptitle(
            f"{bl}: element occurrences, populations and 50 m demes. "
            "Zoom panels: one per population;\ndeme outlines dashed + open "
            "circles = 2025, solid + filled circles = 2026; circle size ∝ "
            "fertile plants per slick spot. Polygons = convex hulls "
            "buffered 25 m (element occurrences 50 m).", fontsize=10)
        out = FIGURES / f"step30h_BL_map_{bl}.png"
        fig.savefig(out, dpi=170, bbox_inches="tight")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
        print(f"[step09] Wrote {out.name}")
    pdf.close()
    print("[step09] Wrote step30h_BL_maps.pdf (all BLs)")


if __name__ == "__main__":
    main()
