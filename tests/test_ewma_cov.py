"""EWMA (RiskMetrics) covariance and correlation matrices."""

import math
import random

import pytest

from quantforge import ewma_covariance_matrix, ewma_correlation_matrix


def _returns(seed=1):
    rng = random.Random(seed)
    out = []
    for _ in range(500):
        z = rng.gauss(0, 1)
        out.append([0.01 * z + 0.005 * rng.gauss(0, 1),
                    0.01 * z + 0.005 * rng.gauss(0, 1),
                    0.02 * rng.gauss(0, 1)])
    return out


def _pd(M):
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


def test_symmetric_and_positive_definite():
    cov = ewma_covariance_matrix(_returns(), 0.94)
    assert all(abs(cov[i][j] - cov[j][i]) < 1e-15 for i in range(3) for j in range(3))
    assert _pd(cov)


def test_correlation_unit_diagonal_and_bounds():
    corr = ewma_correlation_matrix(_returns(), 0.94)
    assert all(abs(corr[i][i] - 1.0) < 1e-9 for i in range(3))
    assert all(-1.0 <= corr[i][j] <= 1.0 for i in range(3) for j in range(3))


def test_shared_factor_positive_independent_zero():
    corr = ewma_correlation_matrix(_returns(), 0.94)
    assert corr[0][1] > 0.5          # assets 0,1 share a factor
    assert abs(corr[0][2]) < 0.2     # asset 2 independent


def test_validation():
    r = _returns()
    with pytest.raises(ValueError):
        ewma_covariance_matrix(r, 1.0)          # lam out of range
    with pytest.raises(ValueError):
        ewma_covariance_matrix([[1, 2]], 0.94)  # < 2 rows
    with pytest.raises(ValueError):
        ewma_covariance_matrix([[1, 2], [3]], 0.94)  # ragged
