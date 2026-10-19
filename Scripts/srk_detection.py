"""srk_detection.py — how many mothers does a deme need? (empirical prior)

One model for the sampling design (Steps 17, 20) and for the sampling
test (Step 21), the same one the SRK predictions use (Step 22):

1. **Deme pool.** A deme of N adults carries 4·N SRK allele copies drawn
   from the species-wide Fg prior (P1, `step26i_L1_carrier_inventory`)
   — drift in a finite deme, so rare Fgs are often absent.
2. **Sampling.** M mothers with s genotyped seedlings each deliver
   A = 4·M + 2·M·s allele draws from that pool (own genotype +
   paternal alleles), so rare Fgs present in the pool may be missed.
3. **Coverage** = Fgs detected / Fgs present, averaged over pools.

`mothers_for_coverage(N)` is the smallest M (≤ N) whose mean coverage
reaches 0.90 with the design 15 seedlings per mother. Given a pool,
the expected number of Fgs detected is Σ_j [1 − (1 − p_j)^A], so the
mean coverage is exact over the simulated pools (no second Monte Carlo
layer). This replaces the earlier uniform-frequency coupon-collector
rule (all Fgs equally common), which under-sizes large demes.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np

PRIOR_TSV = Path("tables") / "step26i_L1_carrier_inventory.tsv"
PLOIDY = 4
PATERNAL_ALLELES_PER_SEEDLING = 2
SEEDLINGS_PER_MOTHER = 15
COVERAGE_TARGET = 0.90
N_POOLS = 2000
RNG_SEED = 2032


@lru_cache(maxsize=1)
def species_prior() -> np.ndarray:
    """P1 species-wide Fg frequencies (same as Step 22)."""
    from step22_predictions_core import build_p1_prior
    return build_p1_prior(PRIOR_TSV)["f_mean"].to_numpy()


@lru_cache(maxsize=None)
def _pool_freqs(N: int) -> tuple[np.ndarray, np.ndarray]:
    """N_POOLS simulated deme pools of N adults: per-pool Fg
    frequencies (N_POOLS, K) and number of Fgs present."""
    f = species_prior()
    rng = np.random.default_rng([RNG_SEED, int(N)])
    counts = rng.multinomial(PLOIDY * int(N), f, size=N_POOLS)
    p = counts / (PLOIDY * int(N))
    return p, (counts > 0).sum(axis=1)


def deme_coverage(N: int, A: int) -> float:
    """Mean Fg coverage of an N-adult deme after A allele draws."""
    if N <= 0:
        return 1.0
    if A <= 0:
        return 0.0
    p, n_present = _pool_freqs(int(N))
    detected = (1.0 - (1.0 - p) ** int(A)).sum(axis=1)
    return float(np.mean(detected / n_present))


def draws(M: int, seedlings_per_mother: int = SEEDLINGS_PER_MOTHER) -> int:
    return M * (PLOIDY + PATERNAL_ALLELES_PER_SEEDLING * seedlings_per_mother)


@lru_cache(maxsize=None)
def mothers_for_coverage(N: int, target: float = COVERAGE_TARGET) -> int:
    """Smallest M ≤ N reaching `target` mean coverage (15 seedlings per
    mother); N when even every adult falls short."""
    N = int(N)
    if N <= 1:
        return max(N, 1)
    if deme_coverage(N, draws(N)) < target:
        return N
    lo, hi = 1, N
    while lo < hi:
        mid = (lo + hi) // 2
        if deme_coverage(N, draws(mid)) >= target:
            hi = mid
        else:
            lo = mid + 1
    return lo
