"""Principal component analysis via Jacobi eigendecomposition."""

import pytest

import math

from quantforge import (
    jacobi_eigen, pca, project, reconstruct_covariance, pca_scenario,
)


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


COV = [[0.04, 0.02, 0.01], [0.02, 0.03, 0.015], [0.01, 0.015, 0.025]]


def test_diagonal_eigenvalues():
    ev, _ = jacobi_eigen([[2, 0, 0], [0, 3, 0], [0, 0, 1]])
    assert ev == sorted(ev, reverse=True)
    assert set(round(x, 6) for x in ev) == {3.0, 2.0, 1.0}


def test_2x2_known_eigenvalues():
    ev, _ = jacobi_eigen([[2, 1], [1, 2]])
    assert ev[0] == pytest.approx(3.0)
    assert ev[1] == pytest.approx(1.0)


def test_eigenvectors_orthonormal():
    _, vec = jacobi_eigen([[2, 1], [1, 2]])
    assert _dot(vec[0], vec[0]) == pytest.approx(1.0)
    assert _dot(vec[0], vec[1]) == pytest.approx(0.0, abs=1e-9)


def test_variances_sum_to_trace():
    res = pca(COV)
    assert sum(res["variances"]) == pytest.approx(0.04 + 0.03 + 0.025)


def test_explained_and_cumulative():
    res = pca(COV)
    assert sum(res["explained"]) == pytest.approx(1.0)
    assert res["cumulative_explained"][-1] == pytest.approx(1.0)
    assert all(res["variances"][i] >= res["variances"][i + 1] for i in range(2))


def test_eigen_equation():
    res = pca(COV)
    v0, l0 = res["loadings"][0], res["variances"][0]
    Av = [sum(COV[i][j] * v0[j] for j in range(3)) for i in range(3)]
    assert all(Av[i] == pytest.approx(l0 * v0[i], abs=1e-9) for i in range(3))


def test_projection_preserves_norm():
    res = pca(COV)
    x = [0.01, -0.02, 0.015]
    scores = project(x, res["loadings"])
    assert _dot(scores, scores) == pytest.approx(_dot(x, x), abs=1e-9)


def test_full_reconstruction_recovers_covariance():
    res = pca(COV)
    full = reconstruct_covariance(res["variances"], res["loadings"])
    assert all(full[i][j] == pytest.approx(COV[i][j], abs=1e-9)
               for i in range(3) for j in range(3))


def test_rank_one_trace_is_top_variance():
    res = pca(COV)
    r1 = reconstruct_covariance(res["variances"], res["loadings"], k=1)
    assert sum(r1[i][i] for i in range(3)) == pytest.approx(res["variances"][0])
    assert all(r1[i][j] == pytest.approx(r1[j][i]) for i in range(3) for j in range(3))


def test_scenario_magnitude_and_scaling():
    res = pca(COV)
    sh = pca_scenario(0, 1.0, res["variances"], res["loadings"])
    mag = math.sqrt(sum(x * x for x in sh))
    assert mag == pytest.approx(math.sqrt(res["variances"][0]))
    sh2 = pca_scenario(0, 2.0, res["variances"], res["loadings"])
    assert all(sh2[i] == pytest.approx(2 * sh[i]) for i in range(3))
    shn = pca_scenario(0, -1.0, res["variances"], res["loadings"])
    assert all(shn[i] == pytest.approx(-sh[i]) for i in range(3))


def test_scenario_validation():
    res = pca(COV)
    with pytest.raises(ValueError):
        pca_scenario(5, 1.0, res["variances"], res["loadings"])


def test_validation():
    with pytest.raises(ValueError):
        jacobi_eigen([[1, 2], [3, 4]])   # asymmetric
