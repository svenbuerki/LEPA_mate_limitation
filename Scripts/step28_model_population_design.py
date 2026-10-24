"""step28_model_population_design.py — deme-level sampling plan for the
model populations (default: D1 = P36 EO25-B small + P34 EO18-7/8 big).

For each population and field season:

1. **Fit demes.** A deme is dropped when it has a single fertile plant
   (no mate within 50 m, so no compatibility to measure) or when fewer
   usable mothers are held than its target (Step 21: mothers selected
   once every sampled seed lot is cleaned < deme target).
2. **Deme selection.** Each deme keeps its size-scaled target (Step 17:
   mothers for 90 % SRK allele detection, `srk_detection`). Demes are
   added greedily by expected new SRK alleles per mother, subject to a
   minimum separation between selected demes, until the expected share of
   the *whole* population's alleles detected (unfit demes included in the
   denominator) reaches the target.
   - Greedy is run from every fit deme as first pick; the plan reaching
     the target with the fewest mothers is kept (`best_plan`).
   - Robustness: demes are chosen on one set of simulated deme pools and
     the stopping rule uses N_EVAL independent pool sets, so the reported
     coverage is not optimistic; the 5th percentile is reported too.
   - Stability: the selection is repeated on N_STABILITY pool sets and the
     share of runs that pick each deme is recorded.
   - Small populations: if all fit demes together need at most
     SMALL_POP_MOTHERS mothers, all are taken (spacing is moot).
   - Adjustment when the plan falls short (decided 2026-10-23): if the
     original plan cannot reach the target, demes dropped only for holding
     too few mothers are re-admitted with every mother they hold, and the
     selection is rerun. Both plans are reported (`original_*` columns).
     In D1 this applies to P34 2025, whose largest deme (158 plants) has
     9 of the 12 mothers its target asks for.
3. **Pick list.** Within each selected deme, mothers are drawn at random
   (reproducible seed per population, year and deme) among those whose
   cleaned seed lot holds ≥ 25 seeds — one per slick spot, then the
   remaining places from the deme's other mothers; lots with fewer seeds,
   then uncleaned lots, are used only when eligible mothers run out. Seed
   count beyond 25 plays no part, so the sample is not biased towards the
   most fecund mothers. Each mother gets 25 seeds to germinate (≈ 15
   seedlings).

Demes are simulated independently, plant by plant, from the species-wide
SRK frequencies (`srk_detection`: tetraploid empirical zygosity, fathers
compatible under sporophytic Class I / II SI); real demes have drifted
further and share alleles with neighbours (Part C), so spatial
separation is a safeguard beyond what the model rewards.

Outputs (no coordinates): tables/step30_D_model_design_demes.tsv,
tables/step30_D_model_design_germplasm_by_deme.csv (lab list),
tables/step30_D_model_design_lab_sheet.tsv (lab sheet keyed on occurrenceID),
tables/step30_D_model_design_summary.tsv,
tables/step30_D_model_design_pick_list.tsv,
figures/step30_D_model_design_curves.{png,pdf}
"""
from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from project_paths import LEPA_DB
from srk_bl_constants import BL_COLORS
from srk_detection import SEEDLINGS_PER_MOTHER, deme_detection, popcount
from srk_units import CROSSWALK_TSV, load_event_units

TABLES = Path("tables")
FIGURES = Path("figures")

# ---- design settings (Sven, 2026-10-22) ----
MODEL_PAIRS = {"D1": {"small": 36, "big": 34},
               "D2": {"small": 21, "big": 19},
               "D3": {"small": 6, "big": 16}}
YEARS = (2025, 2026)
COVERAGE_TARGET = 0.90
MIN_SEPARATION_M = 100.0
MIN_PLANTS_PER_DEME = 2
SMALL_POP_MOTHERS = 10
SEEDS_TO_GERMINATE = 25
N_POOLS = 1000              # simulated populations used to choose demes
N_POOLS_EVAL = 250          # per independent evaluation set
N_EVAL = 20
N_STABILITY = 20
N_RANDOM_BASELINE = 100
RNG_SEED = 2035


def deme_pools(sizes: np.ndarray, M: np.ndarray, rng: np.random.Generator,
               n_rep: int = N_POOLS) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """One simulated population: per deme, the Fgs its adults hold and the
    Fgs its M design mothers (15 seedlings each) detect — (n_rep,)
    uint64 masks each."""
    held, det = [], []
    for n, m in zip(sizes, M):
        p, dt = deme_detection(int(n), int(m), n_rep, rng)
        held.append(p); det.append(dt)
    return held, det


def coverage(sel: list[int], pools) -> np.ndarray:
    """Per-simulation share of the population's alleles (union over ALL
    demes) detected by the selected demes at their design effort."""
    held, det = pools
    n = len(held[0])
    if not sel:
        return np.zeros(n)
    present = np.bitwise_or.reduce(np.stack(held), axis=0)
    found = np.bitwise_or.reduce(np.stack([det[i] for i in sel]), axis=0)
    return popcount(found) / popcount(present)


def deme_distances(year: int, pop: int, deme_ids: list[int]) -> np.ndarray:
    """Shortest slick-spot-to-slick-spot distance (m) between demes."""
    cw = pd.read_csv(CROSSWALK_TSV, sep="\t", usecols=["event_year", "eventID", "lat", "lon"])
    u = load_event_units(year)
    ev = u[u["populationID"] == pop].merge(cw, on=["event_year", "eventID"])
    pts = {d: np.radians(ev.loc[ev["demeID"] == d, ["lat", "lon"]].to_numpy()) for d in deme_ids}
    n = len(deme_ids)
    D = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            p, q = pts[deme_ids[a]], pts[deme_ids[b]]
            dlat = q[None, :, 0] - p[:, None, 0]
            dlon = q[None, :, 1] - p[:, None, 1]
            h = (np.sin(dlat / 2) ** 2 + np.cos(p[:, None, 0]) * np.cos(q[None, :, 0])
                 * np.sin(dlon / 2) ** 2)
            D[a, b] = D[b, a] = (2 * 6_371_008.8 * np.arcsin(np.sqrt(h))).min()
    return D


def greedy(fit: list[int], pools_sel, pools_eval, M, D, min_sep, start=None):
    """Add demes by gain per mother (on pools_sel) until the mean coverage
    on the independent pools_eval reaches the target."""
    sel = [] if start is None else [start]
    while True:
        ev = np.concatenate([coverage(sel, pe) for pe in pools_eval]) if sel else np.zeros(1)
        if ev.mean() >= COVERAGE_TARGET:
            break
        cand = [i for i in fit if i not in sel and all(D[i, j] >= min_sep for j in sel)]
        if not cand:
            break
        c0 = coverage(sel, pools_sel).mean()
        sel.append(max(cand, key=lambda i: (coverage(sel + [i], pools_sel).mean() - c0) / M[i]))
    return sel


def best_plan(fit: list[int], pools_sel, pools_eval, M, D, min_sep):
    """Greedy is short-sighted: its first cheap picks can block, through
    the spacing rule, a deme the population needs. Run it from every fit
    deme as first pick and keep the plan that reaches the target with the
    fewest mothers (ties: higher coverage); if none reaches it, the plan
    with the highest coverage."""
    best, key_best = None, None
    for s in [None] + list(fit):
        sel = greedy(fit, pools_sel, pools_eval, M, D, min_sep, start=s)
        cov = float(np.mean([coverage(sel, pe).mean() for pe in pools_eval]))
        ok = cov >= COVERAGE_TARGET
        key = (not ok, int(M[sel].sum()) if ok else -cov, -cov)
        if key_best is None or key < key_best:
            best, key_best = sel, key
    return best


def load_mothers() -> pd.DataFrame:
    con = sqlite3.connect(f"file:{LEPA_DB}?mode=ro", uri=True)
    df = pd.read_sql("""
        SELECT o.eventID, o.occurrenceID,
               MAX(g.germplasmID) AS germplasmID,
               COALESCE(MAX(g.germplasmQuantityEstimate), 0) AS seeds
        FROM Occurrences o JOIN Taxonomy t ON o.taxonID = t.taxonID
        LEFT JOIN Germplasm g ON g.occurrenceID = o.occurrenceID
        WHERE t.genus = 'Lepidium' AND t.specificEpithet = 'papilliferum'
          AND (o.provenance IS NULL OR o.provenance = 'in situ')
          AND o.eventID IS NOT NULL
        GROUP BY o.occurrenceID""", con)
    con.close()
    df["cleaned"] = df["germplasmID"].notna()
    germ = np.minimum(df["seeds"].astype(int), SEEDS_TO_GERMINATE)
    df["seedlings"] = np.where(df["cleaned"],
                               np.minimum(np.rint(germ * 0.60).astype(int), SEEDLINGS_PER_MOTHER),
                               SEEDLINGS_PER_MOTHER)
    return df[df["seedlings"] > 0]


def pick_mothers(m: pd.DataFrame, target: int,
                 rng: np.random.Generator) -> pd.DataFrame:
    """Random mothers within one deme, not ranked by seed yield.

    Tiers: (1) eligible — cleaned seed lot with ≥ SEEDS_TO_GERMINATE seeds
    (enough for the full 15 seedlings); (2) cleaned lot with fewer seeds
    (germinated in full); (3) seed not cleaned yet. One mother per slick
    spot is drawn at random from the best tier available in that slick
    spot; the remaining places are drawn at random from the deme's other
    mothers, best tier first. Seed count beyond the 25-seed threshold plays
    no part, so the sample is not biased towards the most fecund mothers.
    """
    m = m.copy()
    m["tier"] = np.select([m["cleaned"] & (m["seeds"] >= SEEDS_TO_GERMINATE), m["cleaned"]],
                          [1, 2], 3)
    m["draw"] = rng.random(len(m))
    m = m.sort_values(["tier", "draw"])
    floor = m.groupby("eventID", sort=False).head(1).assign(reason="one per slick spot (random)")
    rest = (m.drop(floor.index).head(max(target - len(floor), 0))
             .assign(reason="random among the deme's other mothers"))
    out = pd.concat([floor, rest])
    tier_note = {1: "", 2: "; lot < 25 seeds", 3: "; seed not cleaned yet"}
    out["reason"] = out["reason"] + out["tier"].map(tier_note)
    return out.drop(columns=["tier", "draw"])


def _label(s):
    return f"P{s['pop']} {s['year']}"


def plot_simulations(sims: list[dict], out: Path) -> None:
    """Four panels validating the design with its own simulations."""
    from srk_detection import deme_coverage, mothers_for_coverage
    fig = plt.figure(figsize=(15, 16))
    gs = fig.add_gridspec(3, 2, hspace=0.55, wspace=0.22)

    # A — Step 1: detection vs mothers per deme size
    ax = fig.add_subplot(gs[0, 0])
    Ms = np.arange(1, 31)
    cmap = plt.get_cmap("viridis")
    for k, n in enumerate([5, 10, 20, 50, 100, 200, 400]):
        cov = [deme_coverage(n, m) for m in Ms if m <= n]
        col = cmap(k / 6)
        ax.plot(Ms[:len(cov)], cov, color=col, lw=1.6, label=f"{n} plants")
        mt = mothers_for_coverage(n)
        ax.plot(mt, deme_coverage(n, mt), "o", color=col, mec="black", ms=6)
    ax.axhline(COVERAGE_TARGET, ls=":", c="#444")
    ax.set(xlabel="Mothers sampled in the deme (15 seedlings each)",
           ylabel="Share of the deme's SRK alleles detected", ylim=(0, 1.02))
    ax.set_title("A. Step 1 — mothers needed per deme (dot = target)", loc="left",
                 weight="bold", fontsize=11)
    ax.legend(title="Deme size", fontsize=8, title_fontsize=8, loc="lower right")

    # B — alleles held vs detected in the selected demes
    ax = fig.add_subplot(gs[0, 1])
    x, labels, xpos = 0, [], []
    for s in sims:
        for dd in s["per_deme"]:
            h, dt = dd["held"], dd["detected"]
            ax.errorbar(x - 0.15, h.mean(), yerr=[[h.mean() - np.quantile(h, .025)],
                        [np.quantile(h, .975) - h.mean()]], fmt="s", color="#888", ms=4, capsize=2)
            ax.errorbar(x + 0.15, dt.mean(), yerr=[[dt.mean() - np.quantile(dt, .025)],
                        [np.quantile(dt, .975) - dt.mean()]], fmt="o",
                        color="#E41A1C", ms=4, capsize=2)
            labels.append(f"{_label(s)} d{dd['demeID']} ({dd['plants']} pl., {dd['mothers']} m.)")
            xpos.append(x); x += 1
        x += 0.6
    ax.set_xticks(xpos); ax.set_xticklabels(labels, rotation=70, ha="right", fontsize=7)
    ax.set_ylabel("SRK alleles per deme")
    ax.plot([], [], "s", color="#888", label="held by the deme")
    ax.plot([], [], "o", color="#E41A1C", label="detected by its mothers' seedlings")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title("B. Selected demes — alleles held vs detected (95 % interval)", loc="left",
                 weight="bold", fontsize=11)

    # C — population coverage of the plan vs every fit deme
    ax = fig.add_subplot(gs[1, 0])
    cols = {"selected": "#E41A1C", "all_fit": "#377EB8"}
    ticks, tlabs = [], []
    for k, s in enumerate(sims):
        for j, key in enumerate(["selected", "all_fit"]):
            parts = ax.violinplot([s[key]], positions=[k * 2.4 + j * 0.9], widths=0.8,
                                  showmedians=True, showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(cols[key]); b.set_alpha(0.6)
            parts["cmedians"].set_color("black")
        ticks.append(k * 2.4 + 0.45)
        tlabs.append(f"{_label(s)}\nplan {s['mothers']} · all {s['all_fit_mothers']} mothers")
    ax.axhline(COVERAGE_TARGET, ls=":", c="#444")
    ax.set_xticks(ticks); ax.set_xticklabels(tlabs, fontsize=8)
    ax.set_ylabel("Share of the population's SRK alleles detected\n(independent simulations)")
    ax.set_ylim(0.6, 1.01)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=cols["selected"], alpha=0.6, label="plan (selected demes)"),
                       plt.Rectangle((0, 0), 1, 1, color=cols["all_fit"], alpha=0.6, label="every fit deme")],
              fontsize=8, loc="lower right")
    ax.set_title("C. Population coverage — plan vs sampling every fit deme", loc="left",
                 weight="bold", fontsize=11)

    # D — mothers needed: plan vs demes added in random order
    ax = fig.add_subplot(gs[1, 1])
    for k, s in enumerate(sims):
        rm = s["random_mothers"]
        parts = ax.violinplot([rm], positions=[k], widths=0.7, showextrema=False)
        for b in parts["bodies"]:
            b.set_facecolor("#999999"); b.set_alpha(0.6)
        ax.plot([k - 0.35, k + 0.35], [s["mothers"]] * 2, color="#E41A1C", lw=2.5)
        ax.text(k + 0.38, s["mothers"], f"{s['mothers']}", color="#E41A1C", va="center", fontsize=8)
        ax.text(k, rm.max() * 1.03 + 1, f"{100 * np.mean(rm > s['mothers']):.0f} % of\nrandom orders\nneed more",
                ha="center", fontsize=7.5)
    ax.set_xticks(range(len(sims))); ax.set_xticklabels([_label(s) for s in sims], fontsize=8)
    ax.set_ylabel("Mothers needed to reach 90 % of the\npopulation's SRK alleles")
    ax.set_ylim(0, max(s["random_mothers"].max() for s in sims) * 1.25 + 5)
    ax.legend(handles=[plt.Line2D([], [], color="#E41A1C", lw=2.5, label="plan (Step 2 greedy, ≥ 100 m)"),
                       plt.Rectangle((0, 0), 1, 1, color="#999999", alpha=0.6,
                                     label=f"demes added in random order ({N_RANDOM_BASELINE} runs)")],
              fontsize=8, loc="upper left")
    ax.set_title("D. Step 2 validated — mothers needed: plan vs random deme order", loc="left",
                 weight="bold", fontsize=11)

    # E — stability of the deme selection
    ax = fig.add_subplot(gs[2, :])
    x, xt, xl = 0, [], []
    for s in sims:
        if not s["rule"].startswith("greedy"):
            continue
        for dm, fq in sorted(s["freq"].items()):
            ax.bar(x, fq, color="#E41A1C" if dm in s["sel_ids"] else "#BBBBBB", width=0.8)
            xt.append(x); xl.append(f"{_label(s)} d{dm}"); x += 1
        x += 1
    ax.set_xticks(xt); ax.set_xticklabels(xl, rotation=70, ha="right", fontsize=7)
    ax.set_ylabel("Share of reruns selecting the deme")
    ax.set_ylim(0, 1.05)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color="#E41A1C", label="in the plan"),
                       plt.Rectangle((0, 0), 1, 1, color="#BBBBBB", label="not in the plan")], fontsize=8)
    ax.set_title(f"E. Stability — {N_STABILITY} reruns on fresh simulations (fit demes)", loc="left",
                 weight="bold", fontsize=11)
    for a in fig.axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Model-population sampling design — simulation results", fontsize=13, weight="bold", y=0.93)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_schematic_map(pop: int, demes: pd.DataFrame, out: Path) -> None:
    """Deme map per year in metres from an arbitrary origin: no basemap,
    no latitude/longitude, only a scale bar. Selected demes coloured
    (label: deme · mothers), fit not selected grey, unfit dashed."""
    from scipy.spatial import ConvexHull
    cw = pd.read_csv(CROSSWALK_TSV, sep="\t", usecols=["event_year", "eventID", "lat", "lon"])
    u = load_event_units()
    ev = u[u["populationID"] == pop].merge(cw, on=["event_year", "eventID"])
    lat0 = ev["lat"].mean()
    ev["x"] = (ev["lon"] - ev["lon"].mean()) * 111_320 * np.cos(np.radians(lat0))
    ev["y"] = (ev["lat"] - ev["lat"].mean()) * 110_540
    ev["x"] -= ev["x"].min(); ev["y"] -= ev["y"].min()   # arbitrary origin
    col = BL_COLORS.get(ev["BL"].iat[0], "#555")
    from matplotlib.colors import to_rgb
    label_col = tuple(0.55 * c for c in to_rgb(col))   # darker shade of the fill
    years = sorted(ev["event_year"].unique())
    span = max(np.ptp(ev["x"]), np.ptp(ev["y"]), 200) + 120
    fig, axes = plt.subplots(1, len(years), figsize=(7.5 * len(years), 7.5))
    axes = np.atleast_1d(axes)
    a = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    grey_labels = []          # placed last, away from the sampled labels
    for ax, yr in zip(axes, years):
        e = ev[ev["event_year"] == yr]
        dy = demes[(demes["populationID"] == pop) & (demes["year"] == yr)].set_index("demeID")
        n_sel = int(dy["selected"].sum()); m_sel = int(dy.loc[dy["selected"], "mothers_used"].sum())
        for dm, g in e.groupby("demeID"):
            pts = np.column_stack([(g["x"].to_numpy()[:, None] + 25 * np.cos(a)).ravel(),
                                   (g["y"].to_numpy()[:, None] + 25 * np.sin(a)).ravel()])
            ring = pts[ConvexHull(pts).vertices]; ring = np.vstack([ring, ring[:1]])
            r = dy.loc[dm]
            if r["selected"]:
                ax.fill(ring[:, 0], ring[:, 1], color=col, alpha=0.55, zorder=2)
                ax.plot(ring[:, 0], ring[:, 1], color="black", lw=0.6, zorder=3)
                ax.annotate(f"d{dm} · {int(r['mothers_used'])}{'*' if r['admitted_with_fewer_mothers'] else ''}", (ring[:, 0].max(), ring[:, 1].max()),
                            fontsize=8, weight="bold", color=label_col, zorder=5,
                            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=label_col, lw=0.6))
            elif r["fit"]:
                ax.fill(ring[:, 0], ring[:, 1], color="#BBBBBB", alpha=0.7, zorder=2)
                ax.plot(ring[:, 0], ring[:, 1], color="#555", lw=0.5, zorder=3)
                grey_labels.append((ax, f"d{dm}", ring))
            else:
                ax.plot(ring[:, 0], ring[:, 1], color="#555", lw=0.8, ls="--", zorder=3)
            ax.scatter(g["x"], g["y"], s=6 + 4 * np.sqrt(g["n_fertile"]), color="black",
                       zorder=4, linewidths=0)
        ax.set_xlim(-60, span); ax.set_ylim(-60, span)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("#999")
        bar = 100 if span < 1500 else 500
        ax.plot([0, bar], [-35, -35], color="black", lw=3)
        ax.text(bar / 2, -25, f"{bar} m", ha="center", va="bottom", fontsize=8)
        ax.set_title(f"{yr} — {n_sel} demes sampled, {m_sel} mothers", loc="left", weight="bold")
    handles = [plt.Rectangle((0, 0), 1, 1, color=col, alpha=0.55, label="sampled (deme · mothers)"),
               plt.Rectangle((0, 0), 1, 1, color="#BBBBBB", label="not sampled (deme)"),
               plt.Line2D([], [], color="#555", ls="--", label="single plant, left out"),
               plt.Line2D([], [], ls="", label="* fewer mothers held than target (all used)"),
               plt.Line2D([], [], marker="o", ls="", color="black", ms=4,
                          label="slick spot (size ∝ fertile plants)")]
    axes[-1].legend(handles=handles, fontsize=8, loc="upper right")
    fig.suptitle(f"P{pop} — demes per year (positions in metres from an arbitrary origin; "
                 "no geographic coordinates)", fontsize=12, weight="bold", color=col, y=1.03)
    fig.tight_layout()
    # grey labels: first corner of the deme (bottom-left, top-left,
    # bottom-right, top-right) whose label overlaps no label drawn so far
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    taken = {ax: [t_.get_window_extent(rend).expanded(1.05, 1.15) for t_ in ax.texts]
             for ax in axes}
    for ax, lab, ring in grey_labels:
        x0, x1 = ring[:, 0].min(), ring[:, 0].max()
        y0, y1 = ring[:, 1].min(), ring[:, 1].max()
        cands = [((x0, y0), "right", "top"), ((x0, y1), "right", "bottom"),
                 ((x1, y0), "left", "top"), ((x1, y1), "left", "bottom")]
        for (xy, ha, va) in cands:
            txt = ax.annotate(lab, xy, fontsize=7, color="#888888", style="italic",
                              ha=ha, va=va, zorder=5)
            bb = txt.get_window_extent(rend)
            if not any(bb.overlaps(o) for o in taken[ax]):
                break
            txt.remove()
        else:
            txt = ax.annotate(lab, cands[0][0], fontsize=7, color="#888888", style="italic",
                              ha="right", va="top", zorder=5)
            bb = txt.get_window_extent(rend)
        taken[ax].append(bb.expanded(1.05, 1.15))
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", default="D1", help="comma-separated model pairs (D1, D2, D3)")
    args = ap.parse_args()
    pairs = [p.strip() for p in args.pairs.split(",") if p.strip()]
    mothers = load_mothers()
    deme_rows, summary_rows, pick_rows, curves, sims = [], [], [], [], []
    for pair in pairs:
        for role, pop in MODEL_PAIRS[pair].items():
            for year in YEARS:
                dt = pd.read_csv(TABLES / f"step30_A_sampling_adequacy_per_deme_{year}.tsv", sep="\t")
                d = dt[dt["populationID"] == pop].reset_index(drop=True)
                if d.empty:
                    continue
                many = d["deme_N_fertile"] >= MIN_PLANTS_PER_DEME
                d["fit"] = many & (d["M_selected_capacity"] >= d["deme_M_target"])
                short = many & (d["M_selected_capacity"] >= 1) & ~d["fit"]
                d["unfit_reason"] = np.select(
                    [~many, d["M_selected_capacity"] < d["deme_M_target"]],
                    ["single plant", "too few mothers held"], "")
                d["mothers_used"] = np.minimum(d["deme_M_target"],
                                               d["M_selected_capacity"]).clip(lower=1)
                M = d["mothers_used"].to_numpy(int)
                sizes = d["deme_N_fertile"].to_numpy(int)
                fit = [i for i in range(len(d)) if d.loc[i, "fit"]]
                D = deme_distances(year, pop, d["demeID"].astype(int).tolist())
                rng = np.random.default_rng([RNG_SEED, pop, year])
                pools_sel = deme_pools(sizes, M, rng)
                pools_eval = [deme_pools(sizes, M, rng, N_POOLS_EVAL) for _ in range(N_EVAL)]

                def plan(cands):
                    if M[cands].sum() <= SMALL_POP_MOTHERS:
                        return list(cands), "all fit demes (small population)"
                    return (best_plan(cands, pools_sel, pools_eval, M, D, MIN_SEPARATION_M),
                            f"greedy (best first deme), ≥ {MIN_SEPARATION_M:.0f} m apart")
                sel, rule = plan(fit)
                orig = dict(demes=len(sel), mothers=int(M[sel].sum()),
                            coverage=float(np.concatenate([coverage(sel, pe) for pe in pools_eval]).mean()))
                d["admitted_short"] = False
                if orig["coverage"] < COVERAGE_TARGET and short.any():
                    d.loc[short, "admitted_short"] = True
                    d.loc[short, "fit"] = True
                    d.loc[short, "unfit_reason"] = ""
                    fit = [i for i in range(len(d)) if d.loc[i, "fit"]]
                    sel, rule = plan(fit)
                    rule += " + short-handed demes admitted"
                cov_eval = np.concatenate([coverage(sel, pe) for pe in pools_eval])
                cov_all = np.concatenate([coverage(fit, pe) for pe in pools_eval])
                # Validation: demes added in random order (same 90 % stopping
                # rule as the plan, no spacing) — how many mothers does that take?
                m_sel = int(M[sel].sum())
                rr = np.random.default_rng([RNG_SEED, pop, year, 99])
                rand_mothers = []
                for b in range(N_RANDOM_BASELINE):
                    pick, mm = [], 0
                    for i in rr.permutation(fit):
                        pick.append(int(i)); mm += int(M[i])
                        if np.mean([coverage(pick, pe).mean() for pe in pools_eval]) >= COVERAGE_TARGET:
                            break
                    rand_mothers.append(mm)
                rand_mothers = np.array(rand_mothers)
                # per selected deme: alleles held vs detected (first evaluation set)
                per_deme = []
                for i in sel:
                    held = popcount(pools_eval[0][0][i])
                    det = popcount(pools_eval[0][1][i])
                    per_deme.append(dict(demeID=int(d.loc[i, "demeID"]), plants=int(sizes[i]),
                                         mothers=int(M[i]), held=held, detected=det))
                freq = np.zeros(len(d))
                if rule.startswith("greedy"):
                    for s in range(N_STABILITY):
                        r2 = np.random.default_rng([RNG_SEED, pop, year, s + 1])
                        ps = deme_pools(sizes, M, r2, N_POOLS_EVAL * 2)
                        for i in best_plan(fit, ps, [deme_pools(sizes, M, r2, N_POOLS_EVAL) for _ in range(3)],
                                           M, D, MIN_SEPARATION_M):
                            freq[i] += 1
                    freq /= N_STABILITY
                else:
                    freq[sel] = 1.0
                order = {i: k + 1 for k, i in enumerate(sel)}
                for i in range(len(d)):
                    r = d.loc[i]
                    deme_rows.append(dict(pair=pair, role=role, populationID=pop, year=year,
                                          demeID=int(r["demeID"]), fertile_plants=int(r["deme_N_fertile"]),
                                          mother_target=int(r["deme_M_target"]),
                                          mothers_used=int(r["mothers_used"]),
                                          admitted_with_fewer_mothers=bool(r["admitted_short"]),
                                          mothers_held=int(r["n_occurrences"]),
                                          mothers_cleaned=int(r["n_cleaned"]),
                                          fit=bool(r["fit"]), unfit_reason=r["unfit_reason"],
                                          selected=i in order, selection_order=order.get(i),
                                          selected_in_share_of_runs=round(freq[i], 2)))
                cum, cm = [], 0
                for k, i in enumerate(sel):
                    cm += M[i]
                    cum.append((cm, float(np.concatenate([coverage(sel[:k + 1], pe)
                                                          for pe in pools_eval]).mean())))
                curves.append(dict(pair=pair, role=role, pop=pop, year=year, cum=cum,
                                   all_fit=(int(M[fit].sum()), float(cov_all.mean()))))
                sims.append(dict(pair=pair, role=role, pop=pop, year=year, mothers=m_sel,
                                 selected=cov_eval, all_fit=cov_all,
                                 random_mothers=rand_mothers,
                                 all_fit_mothers=int(M[fit].sum()), per_deme=per_deme,
                                 freq={int(d.loc[i, "demeID"]): freq[i] for i in fit},
                                 sel_ids={int(d.loc[i, "demeID"]) for i in sel},
                                 rule=rule))
                summary_rows.append(dict(
                    pair=pair, role=role, populationID=pop, year=year, rule=rule,
                    original_selected_demes=orig["demes"], original_mothers=orig["mothers"],
                    original_coverage_mean=round(orig["coverage"], 3),
                    demes=len(d), unfit_demes=int((~d["fit"]).sum()), fit_demes=len(fit),
                    mothers_all_fit=int(M[fit].sum()),
                    coverage_all_fit=round(float(cov_all.mean()), 3),
                    selected_demes=len(sel), mothers_selected=int(M[sel].sum()),
                    coverage_selected_mean=round(float(cov_eval.mean()), 3),
                    coverage_selected_p05=round(float(np.quantile(cov_eval, 0.05)), 3),
                    random_order_mothers_median=int(np.median(rand_mothers)),
                    random_order_mothers_p05_p95=f"{int(np.quantile(rand_mothers, .05))}-{int(np.quantile(rand_mothers, .95))}",
                    random_order_share_needing_more=round(float(np.mean(rand_mothers > m_sel)), 2),
                    seeds_to_germinate=int(M[sel].sum()) * SEEDS_TO_GERMINATE,
                    seedlings_expected=int(M[sel].sum()) * SEEDLINGS_PER_MOTHER))
                # pick list
                u = load_event_units(year)
                mo = mothers.merge(u[["eventID", "populationID", "demeID"]], on="eventID")
                for i in sel:
                    dm = mo[(mo["populationID"] == pop) & (mo["demeID"] == d.loc[i, "demeID"])]
                    rng_pick = np.random.default_rng([RNG_SEED, pop, year,
                                                      int(d.loc[i, "demeID"]), 7])
                    for _, r in pick_mothers(dm, int(M[i]), rng_pick).iterrows():
                        pick_rows.append(dict(
                            pair=pair, role=role, populationID=pop, year=year,
                            demeID=int(d.loc[i, "demeID"]), eventID=int(r["eventID"]),
                            occurrenceID=int(r["occurrenceID"]),
                            germplasmID=(int(r["germplasmID"]) if r["cleaned"] else None),
                            status="germinate" if r["cleaned"] else "clean seed first",
                            seeds_available=int(r["seeds"]) if r["cleaned"] else None,
                            seeds_to_germinate=(min(int(r["seeds"]), SEEDS_TO_GERMINATE)
                                                if r["cleaned"] else SEEDS_TO_GERMINATE),
                            seedlings_expected=int(r["seedlings"]), reason=r["reason"]))
    demes = pd.DataFrame(deme_rows)
    summary = pd.DataFrame(summary_rows)
    picks = pd.DataFrame(pick_rows)
    picks["germplasmID"] = picks["germplasmID"].astype("Int64")
    picks["seeds_available"] = picks["seeds_available"].astype("Int64")
    demes.to_csv(TABLES / "step30_D_model_design_demes.tsv", sep="\t", index=False)
    summary.to_csv(TABLES / "step30_D_model_design_summary.tsv", sep="\t", index=False)
    picks.to_csv(TABLES / "step30_D_model_design_pick_list.tsv", sep="\t", index=False)
    # Lab list: the germplasmIDs to pull, per deme (CSV for the seed lab)
    # small population first, then big — the order used in every document
    role_order = {"small": 0, "big": 1}
    lab = (picks.assign(_r=picks["role"].map(role_order))
                .sort_values(["pair", "_r", "year", "demeID", "status", "germplasmID"])
                .drop(columns="_r"))
    lab = lab.assign(population=lambda x: "P" + x["populationID"].astype(str))[
        ["pair", "role", "population", "year", "demeID", "germplasmID", "occurrenceID",
         "status", "seeds_available", "seeds_to_germinate", "seedlings_expected", "reason"]]
    lab.to_csv(TABLES / "step30_D_model_design_germplasm_by_deme.csv", index=False)
    # Lab sheet (TSV) keyed on occurrenceID — germplasmIDs are still missing
    # for seed lots not yet cleaned — with BL, population, EO and deme.
    u_all = load_event_units()
    ev_info = u_all.drop_duplicates("eventID")[["eventID", "BL", "locationCode"]]
    sheet = picks.merge(ev_info, on="eventID", how="left")
    sheet["EO"] = sheet["locationCode"].str.extract(r"^(EO\d+)")[0].fillna(sheet["locationCode"])
    sheet = sheet.assign(population="P" + sheet["populationID"].astype(str))[
        ["occurrenceID", "germplasmID", "BL", "population", "EO", "locationCode", "year",
         "demeID", "eventID", "pair", "role", "status", "seeds_available",
         "seeds_to_germinate", "seedlings_expected", "reason"]]
    sheet = sheet.rename(columns={"eventID": "eventID_slick_spot", "reason": "selected_because"})
    sheet = (sheet.assign(_r=sheet["role"].map(role_order))
                  .sort_values(["BL", "_r", "year", "demeID", "eventID_slick_spot", "occurrenceID"])
                  .drop(columns="_r"))
    sheet.to_csv(TABLES / "step30_D_model_design_lab_sheet.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print(f"[step28] pick list: {len(picks)} mothers "
          f"({(picks['status'] == 'germinate').sum()} ready, "
          f"{(picks['status'] != 'germinate').sum()} to clean first)")

    fig, axes = plt.subplots(1, len(curves), figsize=(3.6 * len(curves), 3.8), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, c in zip(axes, curves):
        uy = load_event_units(c["year"])
        col = BL_COLORS.get(uy.loc[uy["populationID"] == c["pop"], "BL"].iat[0], "#555")
        if c["cum"]:
            x, y = zip(*([(0, 0.0)] + c["cum"]))
            ax.plot(x, y, "-o", color=col, ms=4)
        ax.plot(*c["all_fit"], "s", color="#444", ms=6, mfc="white")
        ax.annotate("all fit demes", c["all_fit"], xytext=(4, -12),
                    textcoords="offset points", fontsize=7)
        ax.axhline(COVERAGE_TARGET, ls=":", c="#444", lw=0.8)
        ax.set_ylim(0, 1.02)
        ax.set_title(f"{c['pair']} {c['role']} — P{c['pop']}, {c['year']}", fontsize=9)
        ax.set_xlabel("Mothers sampled (cumulative)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Share of the population's SRK alleles\ndetected (independent simulations)")
    fig.suptitle("Model-population sampling: demes added one at a time (best gain per "
                 "mother, ≥ 100 m apart) until 90 % of the population's SRK alleles "
                 "are expected to be detected", fontsize=9)
    fig.tight_layout()
    FIGURES.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"step30_D_model_design_curves.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    plot_simulations(sims, FIGURES / "step30_D_model_design_simulations")
    for pop in sorted({s["pop"] for s in sims}):
        plot_schematic_map(pop, demes, FIGURES / f"step30_D_model_design_demes_schematic_P{pop}")


if __name__ == "__main__":
    main()
