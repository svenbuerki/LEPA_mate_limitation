"""Step 03 (historical tag: Step 30h) — hierarchical geographic
clustering of the populations into Bottleneck Lineages (BLs).

§ A.4.5 (new). Replicates the content of the external
LEPA_EO_spatial_clustering repo (Ward's D2 at the EO level) but now
at the POPULATION level, using the 2 years of pooled 2025+2026 data
that produced the 44 populations in step29a. The purpose is to
re-check whether the Bottleneck-Lineage (BL) grouping used
throughout still holds.

**The number of BLs is fixed at `BL_K = 5`** (decision 2026-10-17), in
continuity with the EO-level Bottleneck Lineages. The silhouette curve
is still computed and its optimum reported as a diagnostic: on the
corrected 2025 + 2026 data it is k = 8 (0.747 vs 0.716 at k = 5), on a
plateau from k = 5 to 9, so the extra splits are sub-groups within
lineages.

Approach
--------
1. Compute per-population centroid (lat, lon) from the
   step29a_population_crosswalk.tsv — a size-weighted mean of event
   coordinates across both years so that the centroid tracks the
   bulk of the above-ground sample.
2. Build the pairwise haversine distance matrix (metres).
3. Ward's D2 linkage on that distance matrix (standard choice for
   geographic hierarchical clustering of this scale, matches the
   external repo).
4. Cut the tree at BL_K = 5 → BL clusters.
5. Silhouette for k = 2 … 10 (diagnostic; optimum reported, not used).
6. Dendrogram with the k = 5 cut highlighted; silhouette curve
   alongside.

Outputs
-------
tables/step30h_population_distances.tsv
    long-form distance matrix (populationA, populationB,
    distance_m, pair label).
tables/step30h_cluster_assignments.tsv
    populationID → BL cluster at BL_K = 5 in `cluster_k_optimal` (the
    column Steps 04 – 05 read; the name is historical), the same in
    `cluster_k5`, and the silhouette optimum in `k_silhouette_optimum`
    with its clusters in `cluster_k_silhouette`.
tables/step30h_silhouette_curve.tsv
    k vs silhouette coefficient; the whole curve for the report.
figures/step30h_dendrogram.png / .pdf
figures/step30h_silhouette_curve.png / .pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score

from step14_seed_sampling_per_mother import haversine_meters

TABLES = Path("tables")
FIGURES = Path("figures")

K_RANGE = range(2, 11)
# Number of Bottleneck Lineages. Fixed (decision 2026-10-17) for
# continuity with the EO-level lineages; the silhouette optimum is only
# reported. See the module docstring.
BL_K = 5


# ---------------------------------------------------------------------------
# 1. Centroids
# ---------------------------------------------------------------------------
def load_population_centroids() -> pd.DataFrame:
    """Per-population centroid = size-weighted mean of event (lat, lon)
    across both years, so the centroid sits where the plants actually
    are rather than at the midpoint of a sparsely-sampled hull."""
    cw = pd.read_csv(TABLES / "step29a_population_crosswalk.tsv",
                      sep="\t", encoding="utf-8-sig")
    # weight by n_fertile (fall back to 1 if NaN or 0)
    cw["w"] = cw["n_fertile"].fillna(0).clip(lower=0)
    cw.loc[cw["w"] == 0, "w"] = 1
    rows = []
    for pop_id, g in cw.groupby("populationID"):
        rows.append({
            "populationID": int(pop_id),
            "lat": float(np.average(g["lat"], weights=g["w"])),
            "lon": float(np.average(g["lon"], weights=g["w"])),
            "n_fert_total": int(g["n_fertile"].sum()),
            "n_events": int(g["eventID"].nunique()),
        })
    return pd.DataFrame(rows).sort_values("populationID").reset_index(drop=True)


# ---------------------------------------------------------------------------
# 2. Distances + linkage
# ---------------------------------------------------------------------------
def pairwise_distance_matrix(cents: pd.DataFrame) -> np.ndarray:
    lat = cents["lat"].to_numpy()
    lon = cents["lon"].to_numpy()
    return haversine_meters(
        lat[:, None], lon[:, None],
        lat[None, :], lon[None, :],
    )


def long_form_distances(cents: pd.DataFrame,
                          D: np.ndarray) -> pd.DataFrame:
    rows = []
    ids = cents["populationID"].to_numpy()
    n = len(ids)
    for i in range(n):
        for j in range(i + 1, n):
            rows.append({
                "populationA": int(ids[i]),
                "populationB": int(ids[j]),
                "distance_m": float(D[i, j]),
            })
    df = pd.DataFrame(rows).sort_values("distance_m").reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# 3. Silhouette for k selection
# ---------------------------------------------------------------------------
def silhouette_curve(D: np.ndarray, Z: np.ndarray,
                       k_range) -> pd.DataFrame:
    rows = []
    for k in k_range:
        labels = fcluster(Z, t=k, criterion="maxclust")
        if len(set(labels)) < 2:
            s = np.nan
        else:
            s = silhouette_score(D, labels, metric="precomputed")
        rows.append({"k": int(k), "silhouette": float(s),
                      "n_clusters_actual": int(len(set(labels)))})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 4. Figures
# ---------------------------------------------------------------------------
def plot_dendrogram(Z: np.ndarray, cents: pd.DataFrame,
                      k_optimal: int,
                      out_png: Path, out_pdf: Path) -> dict:
    """Render the dendrogram and ALSO return the leaf order so the
    downstream Phase III script can renumber populations by dendrogram
    position (the user-selected within-BL sort key). The dendrogram
    `leaves` attribute is a list of positional indices into `cents`
    in the x-axis order the plot shows them."""
    fig, ax = plt.subplots(figsize=(14, 7))
    labels = [f"P{int(r):>2}" for r in cents["populationID"]]
    ddata = dendrogram(
        Z,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=8,
        color_threshold=Z[-k_optimal + 1, 2] if k_optimal > 1 else 0,
        ax=ax,
    )
    # Cut line
    if k_optimal > 1:
        cut = Z[-k_optimal + 1, 2]
        ax.axhline(cut, color="#444", linestyle="--", linewidth=1.2,
                   label=f"k = {k_optimal} cut (Bottleneck Lineages)")
        ax.legend(loc="upper right", fontsize=10)
    ax.set_ylabel("Haversine linkage distance (m, Ward's D2)")
    ax.set_title(
        f"Population hierarchical clustering — Ward's D2 on pairwise "
        f"haversine distances\n({len(cents)} populations, pooled 2025 + 2026 data; "
        f"k = {k_optimal} Bottleneck Lineages)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30h] Wrote {out_png.name} + .pdf")
    return ddata


def plot_silhouette_curve(sil: pd.DataFrame, k_optimal: int,
                            out_png: Path, out_pdf: Path,
                            k_bl: int = BL_K) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(sil["k"], sil["silhouette"], "-o",
            color="#5e3c99", linewidth=2, markersize=7)
    ax.axvline(k_optimal, color="#D55E00", linestyle="--", linewidth=1.5,
               label=f"silhouette optimum k = {k_optimal}")
    ax.axvline(k_bl, color="#0072B2", linestyle=":", linewidth=1.6,
               label=f"k = {k_bl} (Bottleneck Lineages used)")
    ax.set_xlabel("Number of clusters (k)")
    ax.set_ylabel("Mean silhouette coefficient")
    ax.set_title(
        "Silhouette-based k selection for Ward's D2 population "
        f"clustering (n = {int(sil.attrs.get('n', 0)) or '—'} populations)\nOptimal k = {k_optimal} "
        f"(silhouette = {sil.loc[sil['k'] == k_optimal, 'silhouette'].iat[0]:.3f})",
        fontsize=11,
    )
    ax.set_xticks(list(sil["k"]))
    ax.legend(loc="best", fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30h] Wrote {out_png.name} + .pdf")


# ---------------------------------------------------------------------------
# 5. Driver
# ---------------------------------------------------------------------------
def main() -> None:
    cents = load_population_centroids()
    print(f"[step30h] Loaded {len(cents)} population centroids")

    D = pairwise_distance_matrix(cents)
    long_D = long_form_distances(cents, D)
    long_D.to_csv(TABLES / "step30h_population_distances.tsv",
                   sep="\t", index=False)
    print(f"[step30h] Wrote step30h_population_distances.tsv "
          f"(nearest pair {long_D.iloc[0]['distance_m']:.0f} m; "
          f"farthest pair {long_D.iloc[-1]['distance_m'] / 1000:.1f} km)")

    # Linkage expects a condensed distance matrix (upper triangle flat)
    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method="ward")

    sil = silhouette_curve(D, Z, K_RANGE)
    sil.to_csv(TABLES / "step30h_silhouette_curve.tsv",
                sep="\t", index=False)
    sil.attrs["n"] = len(cents)
    k_sil = int(sil.loc[sil["silhouette"].idxmax(), "k"])
    print(f"[step30h] Silhouette curve:")
    print(sil.to_string(index=False))
    print(f"[step30h] silhouette optimum k = {k_sil} "
          f"(silhouette = {sil.loc[sil['k'] == k_sil, 'silhouette'].iat[0]:.3f}); "
          f"BLs use the fixed BL_K = {BL_K} "
          f"(silhouette = {sil.loc[sil['k'] == BL_K, 'silhouette'].iat[0]:.3f})")

    # Cluster assignments: always expose both an explicit
    # `cluster_k_optimal` column (downstream scripts key on it) AND
    # the `cluster_k5` column for continuity with the old EO BL
    # framework. When k_optimal happens to be 5 the two are
    # identical, which is itself the message.
    labels_bl = fcluster(Z, t=BL_K, criterion="maxclust")
    labels_sil = fcluster(Z, t=k_sil, criterion="maxclust")
    assigns = cents.copy()
    assigns["cluster_k_optimal"] = labels_bl      # read by Steps 04 – 05
    assigns["cluster_k5"] = fcluster(Z, t=5, criterion="maxclust")
    assigns["k_optimal"] = BL_K
    assigns["cluster_k_silhouette"] = labels_sil
    assigns["k_silhouette_optimum"] = k_sil
    k_optimal = BL_K

    # Figures first (so we can read the dendrogram leaf order out
    # before writing cluster_assignments.tsv).
    ddata = plot_dendrogram(
        Z, cents, k_optimal,
        FIGURES / "step30h_dendrogram.png",
        FIGURES / "step30h_dendrogram.pdf",
    )
    # Dendrogram leaf order → populationID → position (0-indexed from
    # left). Each population gets a `dendrogram_leaf_pos` column.
    leaf_positions: dict[int, int] = {}
    for pos, leaf_idx in enumerate(ddata["leaves"]):
        pop_id = int(cents.iloc[leaf_idx]["populationID"])
        leaf_positions[pop_id] = pos
    assigns["dendrogram_leaf_pos"] = assigns["populationID"].map(leaf_positions)
    assigns.to_csv(TABLES / "step30h_cluster_assignments.tsv",
                    sep="\t", index=False)
    print(f"[step30h] Wrote step30h_cluster_assignments.tsv")
    plot_silhouette_curve(
        sil, k_sil,
        FIGURES / "step30h_silhouette_curve.png",
        FIGURES / "step30h_silhouette_curve.pdf",
    )

    # Headline: cluster sizes at both k levels
    print()
    print(f"[step30h] Cluster sizes at k = {BL_K} (Bottleneck Lineages):")
    print(pd.Series(labels_bl).value_counts().sort_index().to_string())
    print(f"[step30h] Cluster sizes at k = {k_sil} (silhouette optimum, diagnostic):")
    print(pd.Series(labels_sil).value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
