"""srk_detection.py — how many mothers does a deme need? (empirical prior,
tetraploid zygosity, sporophytic Class I / II SI)

One model for the sampling design (Steps 17, 20, 28) and for the
sampling test (Step 21), on the same biology as the SRK predictions
(Step 22, `srk_si_model`):

1. **Deme.** A deme of N adults is simulated plant by plant. Each plant
   is tetraploid: its number of distinct SRK alleles (1–4) follows the
   empirical LEPA zygosity (`tables/srk_zygosity_empirical.tsv`: 66 %
   carry one allele in four copies) and its alleles are drawn from the
   species-wide Fg prior (P1, `step26i_L1_carrier_inventory`) — drift
   in a finite deme, so rare Fgs are often absent.
2. **Mothers.** Mothers are drawn at random among the plants that have
   at least one compatible mate in the deme (a plant with none sets no
   seed and is never collected). Sampling M mothers is conditioned on
   the deme having ≥ M such plants — we only collect from plants that
   set seed, so simulated demes where that is impossible are discarded.
3. **Seedlings.** Each seedling has a father drawn at random among the
   plants compatible with its mother under sporophytic Class I / II SI
   (`tables/srk_fg_class.tsv`; Class I dominant over Class II, expressed
   sets must be disjoint) and inherits 2 of his 4 allele copies.
4. **Detection.** A mother's family shows her own alleles plus the
   paternal alleles of her seedlings. Coverage = Fgs detected / Fgs
   present in the deme's adults, averaged over simulated demes.

`mothers_for_coverage(N)` is the smallest M (≤ N) whose mean coverage
reaches 0.90 with the design 15 seedlings per mother. Allele sets are
stored as 64-bit masks (32 Fgs), so whole demes are simulated cheaply.
Coverage curves per deme size are cached in `cache/srk_detection/`
(no coordinates; safe to delete — they are rebuilt).

Earlier versions drew each plant's 4 copies independently and paternal
alleles from the whole pool, ignoring zygosity and SI; that credited a
mother with 34 independent allele draws and under-sized the targets.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np

PRIOR_TSV = Path("tables") / "step26i_L1_carrier_inventory.tsv"
CACHE = Path("cache") / "srk_detection"
MODEL_VERSION = "v3-zygosity-si-seedset"
PLOIDY = 4
PATERNAL_ALLELES_PER_SEEDLING = 2
SEEDLINGS_PER_MOTHER = 15
COVERAGE_TARGET = 0.90
N_POOLS = 1000              # simulated demes per curve / call
RNG_SEED = 2032
# Among two-allele plants: 60 AABB vs 58 AAAB-type (56 AAAB + 2 AAB0),
# from the source column of tables/srk_zygosity_empirical.tsv.
P_AABB_GIVEN_TWO = 60 / 118
_MAX_TRIES = 40             # rejection draws per seedling before exact fallback
_PAIRS = np.array([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
_POP8 = np.array([bin(i).count("1") for i in range(256)], dtype=np.int64)


def popcount(x: np.ndarray) -> np.ndarray:
    """Number of set bits of each uint64 in `x`."""
    x = np.ascontiguousarray(x, dtype=np.uint64)
    return _POP8[x.view(np.uint8)].reshape(*x.shape, 8).sum(axis=-1)


@lru_cache(maxsize=1)
def _model() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Species prior f (K,), Class I mask (K,), zygosity probs (4,)."""
    from step22_predictions_core import build_p1_prior
    from srk_si_model import build_class_i_mask, load_class_map, load_zygosity_dist
    prior = build_p1_prior(PRIOR_TSV)
    f = prior["f_mean"].to_numpy(float)
    class_i = build_class_i_mask(prior["Fg"].astype(str).tolist(), load_class_map())
    return f / f.sum(), class_i, load_zygosity_dist()


def species_prior() -> np.ndarray:
    """P1 species-wide Fg frequencies (same as Step 22)."""
    return _model()[0]


def _genotypes(R: int, N: int, rng: np.random.Generator,
               f: np.ndarray | None = None) -> np.ndarray:
    """(R, N, 4) Fg indices of tetraploid plants under the empirical
    zygosity; distinct alleles drawn from the prior `f` (default: the
    species prior) without replacement (Gumbel top-k = successive
    weighted sampling)."""
    f0, _, zyg = _model()
    f = f0 if f is None else np.asarray(f, float) / np.sum(f)
    with np.errstate(divide="ignore"):
        logf = np.log(f)
    k = rng.choice(4, size=(R, N), p=zyg) + 1
    keys = logf[None, None, :] + rng.gumbel(size=(R, N, len(f)))
    top = np.argsort(-keys, axis=-1)[..., :4]                   # (R, N, 4)
    a, b, c, d = (top[..., i] for i in range(4))
    g = np.stack([a, a, a, a], axis=-1)
    two = k == 2
    aabb = rng.random((R, N)) < P_AABB_GIVEN_TWO
    tri_b = rng.random((R, N)) < 0.5                            # ABBB vs AAAB
    g2 = np.where(aabb[..., None], np.stack([a, a, b, b], -1),
                  np.where(tri_b[..., None], np.stack([a, b, b, b], -1),
                           np.stack([a, a, a, b], -1)))
    g = np.where(two[..., None], g2, g)
    dbl = rng.integers(0, 3, size=(R, N))                       # which allele is doubled
    abc = np.stack([a, b, c], -1)
    extra = np.take_along_axis(abc, dbl[..., None], axis=-1)[..., 0]
    g = np.where((k == 3)[..., None], np.stack([a, b, c, extra], -1), g)
    g = np.where((k == 4)[..., None], np.stack([a, b, c, d], -1), g)
    return g


def deme_allele_counts(N: int, rng: np.random.Generator,
                       f: np.ndarray | None = None) -> np.ndarray:
    """Drift in one deme of N plants: (K,) counts of SRK allele copies
    (4 per plant) under the empirical zygosity — the shared drift model
    of the predictions (Steps 02, 20, 22, 27) and the sampling design.
    Replaces drawing 4·N independent copies, which ignores that two
    thirds of plants carry one allele in four copies and so overstates
    the alleles a deme holds."""
    N = int(N)
    K = len(_model()[0]) if f is None else len(f)
    if N <= 0:
        return np.zeros(K, int)
    g = _genotypes(1, N, rng, f)[0]
    return np.bincount(g.ravel(), minlength=K)


def _bits(idx: np.ndarray) -> np.ndarray:
    """OR of 1 << idx over the last axis → uint64 masks."""
    return np.bitwise_or.reduce(np.left_shift(np.uint64(1), idx.astype(np.uint64)), axis=-1)


def simulate_deme(N: int, n_rep: int, rng: np.random.Generator,
                  seedlings: int = SEEDLINGS_PER_MOTHER
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Simulate `n_rep` demes of N adults and sample every eligible
    mother in random order. Returns
      present (R,) uint64 — Fgs carried by the deme's adults;
      cum     (R, N) uint64 — Fgs detected by the first m mothers at
              column m − 1 (columns beyond the eligible mothers repeat
              the last value);
      n_elig  (R,) int — plants able to set seed (possible mothers)."""
    N = int(N)
    _, class_i, _ = _model()
    ci_bits = np.uint64(_bits(np.where(class_i)[0]))
    present = np.zeros(n_rep, np.uint64)
    cum = np.zeros((n_rep, N), np.uint64)
    n_elig = np.zeros(n_rep, int)
    chunk = max(1, int(2e6 // max(1, N * seedlings * _MAX_TRIES)))
    for r0 in range(0, n_rep, chunk):
        R = min(chunk, n_rep - r0)
        g = _genotypes(R, N, rng)
        own = _bits(g)                                          # (R, N)
        expr = np.where(own & ci_bits, own & ci_bits, own)      # Class I dominance
        # eligibility: at least one compatible plant in the deme
        comp_any = np.empty((R, N), bool)
        for r in range(R):                                      # via distinct expressed sets
            u, inv = np.unique(expr[r], return_inverse=True)
            comp_any[r] = ((u[:, None] & u[None, :]) == 0).any(axis=1)[inv]
        # fathers by rejection sampling, exact fallback for the rest
        cand = rng.integers(0, N, size=(R, N, seedlings, _MAX_TRIES))
        rr = np.arange(R)[:, None, None, None]
        ok = (expr[rr, cand] & expr[:, :, None, None]) == 0
        first = ok.argmax(axis=-1)
        father = np.take_along_axis(cand, first[..., None], axis=-1)[..., 0]
        miss = ~ok.any(axis=-1) & comp_any[:, :, None]
        for r, i, s in zip(*np.nonzero(miss)):
            pool = np.nonzero((expr[r] & expr[r, i]) == 0)[0]
            father[r, i, s] = rng.choice(pool)
        pair = _PAIRS[rng.integers(0, 6, size=(R, N, seedlings))]   # (R, N, s, 2)
        fg = g[np.arange(R)[:, None, None, None], father[..., None], pair]
        pat = _bits(fg.reshape(R, N, -1))
        # a plant with no compatible mate in the deme shows only her own
        # alleles (used only when no deme with enough seed-setters exists,
        # e.g. a lone plant whose seed came from pollen outside the deme)
        det = np.where(comp_any, own | pat, own)
        # random order, eligible mothers first
        order = np.argsort(rng.random((R, N)) + (~comp_any), axis=1)
        det = np.take_along_axis(det, order, axis=1)
        cum[r0:r0 + R] = np.bitwise_or.accumulate(det, axis=1)
        present[r0:r0 + R] = np.bitwise_or.reduce(own, axis=1)
        n_elig[r0:r0 + R] = comp_any.sum(axis=1)
    return present, cum, n_elig


def deme_detection(N: int, M: int, n_rep: int, rng: np.random.Generator,
                   seedlings: int = SEEDLINGS_PER_MOTHER) -> tuple[np.ndarray, np.ndarray]:
    """(present, detected) uint64 masks of `n_rep` simulated demes of N
    adults sampled with M random mothers, conditioned on the deme
    holding ≥ M plants that set seed (fewer only if such demes are
    vanishingly rare, then the most seed-setting plants found)."""
    N, M = int(N), min(int(M), int(N))
    keep_p, keep_d, best = [], [], None
    for _ in range(50):
        p, cum, ne = simulate_deme(N, n_rep, rng, seedlings)
        if M <= 0:
            return p, np.zeros(n_rep, np.uint64)
        ok = ne >= M
        keep_p.append(p[ok]); keep_d.append(cum[ok, M - 1])
        best = (p, cum, ne) if best is None else best
        if sum(len(k) for k in keep_p) >= n_rep:
            return (np.concatenate(keep_p)[:n_rep], np.concatenate(keep_d)[:n_rep])
    p, cum, ne = best
    top = ne.max()
    idx = np.nonzero(ne == top)[0]
    idx = rng.choice(idx, size=n_rep)
    return p[idx], cum[idx, max(top, 1) - 1]


@lru_cache(maxsize=None)
def coverage_curve(N: int, seedlings: int = SEEDLINGS_PER_MOTHER) -> np.ndarray:
    """Mean coverage after m = 1..N mothers (index m − 1), each value over
    the simulated demes with ≥ m seed-setting plants (NaN when none);
    disk-cached."""
    N, seedlings = int(N), int(seedlings)
    path = CACHE / f"{MODEL_VERSION}_N{N}_s{seedlings}_R{N_POOLS}_seed{RNG_SEED}.npy"
    if path.exists():
        return np.load(path)
    present, cum, ne = simulate_deme(N, N_POOLS, np.random.default_rng([RNG_SEED, N, seedlings]),
                                     seedlings)
    cov = popcount(cum) / popcount(present)[:, None]
    valid = ne[:, None] >= np.arange(1, N + 1)[None, :]
    with np.errstate(invalid="ignore"):
        curve = np.where(valid.any(axis=0),
                         (cov * valid).sum(axis=0) / valid.sum(axis=0), np.nan)
    CACHE.mkdir(parents=True, exist_ok=True)
    np.save(path, curve)
    return curve


def deme_coverage(N: int, M: int, seedlings: int = SEEDLINGS_PER_MOTHER) -> float:
    """Mean Fg coverage of an N-adult deme sampled with M mothers."""
    if N <= 0:
        return 1.0
    if M <= 0:
        return 0.0
    return float(coverage_curve(int(N), int(seedlings))[min(int(M), int(N)) - 1])


@lru_cache(maxsize=None)
def mothers_for_coverage(N: int, target: float = COVERAGE_TARGET) -> int:
    """Smallest M ≤ N reaching `target` mean coverage (15 seedlings per
    mother); N when even every adult falls short."""
    N = int(N)
    if N <= 1:
        return max(N, 1)
    reach = np.nonzero(coverage_curve(N) >= target)[0]
    return int(reach[0]) + 1 if len(reach) else N
