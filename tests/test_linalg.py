"""Matrix utilities: Cholesky, correlation repair, correlated draws."""

import pytest

from quantforge import (
    cholesky, is_positive_definite, correlated_normals, nearest_correlation,
    basket_option_mc,
)
from quantforge.pca import jacobi_eigen
from quantforge.multiasset import basket_option


A = [[4, 2, 1], [2, 3, 0.5], [1, 0.5, 2]]


def test_cholesky_reconstructs():
    L = cholesky(A)
    n = 3
    recon = [[sum(L[i][k] * L[j][k] for k in range(n)) for j in range(n)]
             for i in range(n)]
    assert all(recon[i][j] == pytest.approx(A[i][j], abs=1e-9)
               for i in range(n) for j in range(n))


def test_cholesky_lower_triangular():
    L = cholesky(A)
    assert all(L[i][j] == 0 for i in range(3) for j in range(i + 1, 3))


def test_positive_definite_checks():
    assert is_positive_definite(A)
    assert not is_positive_definite([[1, 2], [2, 1]])


def test_correlated_normals_identity():
    z = [0.5, -1.0, 2.0]
    I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    assert correlated_normals(z, I) == z


def test_correlated_normals_transform():
    corr = [[1, 0.5], [0.5, 1]]
    out = correlated_normals([1.0, 1.0], corr)
    L = cholesky(corr)
    assert out[0] == pytest.approx(L[0][0])
    assert out[1] == pytest.approx(L[1][0] + L[1][1])


def test_nearest_correlation_leaves_valid_unchanged():
    valid = [[1, 0.3, 0.2], [0.3, 1, 0.4], [0.2, 0.4, 1]]
    r = nearest_correlation(valid)
    assert all(r[i][j] == pytest.approx(valid[i][j], abs=1e-6)
               for i in range(3) for j in range(3))


def test_nearest_correlation_repairs_indefinite():
    bad = [[1, 0.9, -0.9], [0.9, 1, 0.9], [-0.9, 0.9, 1]]
    assert not is_positive_definite(bad)
    fixed = nearest_correlation(bad)
    assert all(fixed[i][i] == pytest.approx(1.0, abs=1e-6) for i in range(3))
    assert all(fixed[i][j] == pytest.approx(fixed[j][i], abs=1e-9)
               for i in range(3) for j in range(3))
    vals, _ = jacobi_eigen(fixed)
    assert all(v >= -1e-8 for v in vals)   # positive semidefinite


def test_basket_mc_deterministic():
    corr = [[1, 0.3], [0.3, 1]]
    a = basket_option_mc([100, 100], [0.5, 0.5], 100, 1, 0.05, [0.2, 0.25], corr,
                         "call", 5000)
    b = basket_option_mc([100, 100], [0.5, 0.5], 100, 1, 0.05, [0.2, 0.25], corr,
                         "call", 5000)
    assert a == b


def test_basket_mc_higher_correlation_raises_call():
    lo = basket_option_mc([100, 100], [0.5, 0.5], 100, 1, 0.05, [0.2, 0.25],
                          [[1, 0.3], [0.3, 1]], "call", 60000)
    hi = basket_option_mc([100, 100], [0.5, 0.5], 100, 1, 0.05, [0.2, 0.25],
                          [[1, 0.8], [0.8, 1]], "call", 60000)
    assert hi > lo


def test_basket_mc_validation():
    with pytest.raises(ValueError):
        basket_option_mc([100, 100], [0.5, 0.5], 100, 1, 0.05, [0.2, 0.25],
                         [[1, 0.3]], "call")


@pytest.mark.slow
def test_basket_mc_matches_analytic():
    corr = [[1, 0.3], [0.3, 1]]
    an = basket_option(([100, 100]), [0.5, 0.5], 100, 1, 0.05, [0.2, 0.25], 0.3)
    mc = basket_option_mc([100, 100], [0.5, 0.5], 100, 1, 0.05, [0.2, 0.25], corr,
                          "call", 200000)
    assert abs(an - mc) / an < 0.03


def test_validation():
    with pytest.raises(ValueError):
        cholesky([[1, 2], [2, 1]])
