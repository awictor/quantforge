import pytest

from quantforge import observability_matrix, is_observable, observer_gain
from quantforge.lqr import is_controllable
from quantforge.eigen_general import eigenvalues_general


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _tr(X):
    return [[X[i][j] for i in range(len(X))] for j in range(len(X[0]))]


def _error_eigs(A, C, L):
    n = len(A)
    LC = [[L[i][0] * C[0][j] for j in range(n)] for i in range(n)]
    Aerr = [[A[i][j] - LC[i][j] for j in range(n)] for i in range(n)]
    return eigenvalues_general(Aerr)


def test_observability_is_dual_controllability():
    A = [[1.0, 1.0], [0.0, 1.0]]
    C = [[1.0, 0.0]]
    assert is_observable(A, C) == is_controllable(_tr(A), _tr(C))
    assert is_observable(A, C)


def test_observer_places_error_poles():
    A = [[1.0, 1.0], [0.0, 1.0]]
    C = [[1.0, 0.0]]
    poles = [0.2, 0.4]
    L = observer_gain(A, C, poles)
    eigs = sorted(e.real for e in _error_eigs(A, C, L))
    for got, want in zip(eigs, poles):
        assert close(got, want)


def test_unobservable_pair():
    assert not is_observable([[2.0, 0.0], [0.0, 3.0]], [[1.0, 0.0]])


def test_scalar():
    L = observer_gain([[2.0]], [[1.0]], [0.3])
    assert close(L[0][0], 1.7)


def test_observability_matrix_shape():
    Om = observability_matrix([[1.0, 1.0], [0.0, 1.0]], [[1.0, 0.0]])
    assert len(Om) == 2 and len(Om[0]) == 2


def test_three_state_observer():
    A = [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-0.5, 0.3, 0.8]]
    C = [[1.0, 0.0, 0.0]]
    poles = [0.1, 0.2, 0.3]
    L = observer_gain(A, C, poles)
    eigs = sorted(e.real for e in _error_eigs(A, C, L))
    for got, want in zip(eigs, poles):
        assert close(got, want, 1e-5)


def test_multi_output_raises():
    with pytest.raises(ValueError):
        observer_gain([[1, 0], [0, 1]], [[1, 0], [0, 1]], [0.5, 0.5])
