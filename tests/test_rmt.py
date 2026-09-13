"""Random-matrix-theory correlation denoising (eigenvalue clipping)."""

import random

import pytest

from quantforge import marchenko_pastur_edge, clip_correlation_eigenvalues
from quantforge.pca import jacobi_eigen


def _factor_correlation(n=10, t=200, seed=1):
    rng = random.Random(seed)
    data = []
    for _ in range(t):
        f = rng.gauss(0, 1)
        data.append([0.5 * f + rng.gauss(0, 1) for _ in range(n)])
    mean = [sum(data[k][i] for k in range(t)) / t for i in range(n)]
    cov = [[sum((data[k][i] - mean[i]) * (data[k][j] - mean[j]) for k in range(t))
            / (t - 1) for j in range(n)] for i in range(n)]
    sd = [cov[i][i] ** 0.5 for i in range(n)]
    return [[cov[i][j] / (sd[i] * sd[j]) for j in range(n)] for i in range(n)]


def test_mp_edge_formula():
    assert abs(marchenko_pastur_edge(10, 100) - (1 + 0.1 ** 0.5) ** 2) < 1e-12


def test_unit_diagonal():
    corr = _factor_correlation()
    d = clip_correlation_eigenvalues(corr, 200)
    assert all(abs(d[i][i] - 1.0) < 1e-9 for i in range(len(d)))


def test_preserves_signal_eigenvalue():
    corr = _factor_correlation()
    d = clip_correlation_eigenvalues(corr, 200)
    ev0, _ = jacobi_eigen(corr)
    ev1, _ = jacobi_eigen(d)
    assert abs(ev1[0] - ev0[0]) < 0.1        # top (signal) eigenvalue kept


def test_collapses_noise_bulk():
    corr = _factor_correlation()
    d = clip_correlation_eigenvalues(corr, 200)
    ev0, _ = jacobi_eigen(corr)
    ev1, _ = jacobi_eigen(d)
    spread0 = ev0[1] - ev0[-1]
    spread1 = ev1[1] - ev1[-1]
    assert spread1 < spread0          # noise eigenvalues brought together


def test_trace_preserved():
    corr = _factor_correlation()
    d = clip_correlation_eigenvalues(corr, 200)
    n = len(corr)
    assert abs(sum(d[i][i] for i in range(n)) - n) < 1e-6


def test_symmetric():
    corr = _factor_correlation()
    d = clip_correlation_eigenvalues(corr, 200)
    n = len(d)
    assert all(abs(d[i][j] - d[j][i]) < 1e-9 for i in range(n) for j in range(n))


def test_validation():
    with pytest.raises(ValueError):
        marchenko_pastur_edge(0, 100)
    with pytest.raises(ValueError):
        clip_correlation_eigenvalues([[1.0, 0.0]], 100)   # not square
    with pytest.raises(ValueError):
        clip_correlation_eigenvalues([[1.0]], 0)
