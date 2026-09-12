"""Ledoit-Wolf shrinkage covariance estimation."""

import math
import random

import pytest

from quantforge import (
    ledoit_wolf_shrinkage, sample_covariance, constant_correlation_target,
)


def _gen(n, seed):
    """Heterogeneous correlation: assets 0,1 share a factor; 2,3 independent."""
    rng = random.Random(seed)
    rows = []
    for _ in range(n):
        f = rng.gauss(0, 1)
        rows.append([
            0.9 * f + 0.2 * rng.gauss(0, 1),
            0.9 * f + 0.2 * rng.gauss(0, 1),
            rng.gauss(0, 1),
            rng.gauss(0, 1),
        ])
    return rows


def _chol_pd(M):
    p = len(M)
    L = [[0.0] * p for _ in range(p)]
    for a in range(p):
        for b in range(a + 1):
            s = M[a][b] - sum(L[a][k] * L[b][k] for k in range(b))
            if a == b:
                if s <= 0:
                    return False
                L[a][b] = math.sqrt(s)
            else:
                L[a][b] = s / L[b][b]
    return True


def test_delta_in_unit_interval():
    _, delta = ledoit_wolf_shrinkage(_gen(50, 1))
    assert 0.0 <= delta <= 1.0


def test_delta_decreases_with_sample_size():
    deltas = [ledoit_wolf_shrinkage(_gen(n, 7))[1] for n in (50, 500, 5000)]
    assert deltas[0] > deltas[1] > deltas[2]


def test_shrinkage_is_convex_combination():
    data = _gen(50, 3)
    sigma, delta = ledoit_wolf_shrinkage(data)
    S = sample_covariance(data)
    F, _ = constant_correlation_target(S)
    for i in range(len(S)):
        for j in range(len(S)):
            expect = delta * F[i][j] + (1 - delta) * S[i][j]
            assert abs(sigma[i][j] - expect) < 1e-12


def test_diagonal_preserved():
    # F_ii = S_ii, so the shrunk diagonal equals the sample variances.
    data = _gen(50, 5)
    sigma, _ = ledoit_wolf_shrinkage(data)
    S = sample_covariance(data)
    for i in range(len(S)):
        assert abs(sigma[i][i] - S[i][i]) < 1e-12


def test_offdiagonal_pulled_toward_target():
    data = _gen(30, 11)
    sigma, _ = ledoit_wolf_shrinkage(data)
    S = sample_covariance(data)
    F, _ = constant_correlation_target(S)
    i, j = 0, 2
    lo, hi = sorted((S[i][j], F[i][j]))
    assert lo - 1e-12 <= sigma[i][j] <= hi + 1e-12


def test_symmetric_and_positive_definite():
    sigma, _ = ledoit_wolf_shrinkage(_gen(20, 13))
    p = len(sigma)
    assert all(abs(sigma[a][b] - sigma[b][a]) < 1e-12
               for a in range(p) for b in range(p))
    assert _chol_pd(sigma)


def test_validation():
    with pytest.raises(ValueError):
        ledoit_wolf_shrinkage([[1.0, 2.0]])  # < 2 observations
    with pytest.raises(ValueError):
        sample_covariance([[1.0, 2.0], [3.0]])  # ragged
