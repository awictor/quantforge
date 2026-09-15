import pytest

from quantforge import ackermann
from quantforge.eigen_general import eigenvalues_general


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _closed_loop_eigs(A, B, K):
    n = len(A)
    BK = [[B[i][0] * K[0][j] for j in range(n)] for i in range(n)]
    Acl = [[A[i][j] - BK[i][j] for j in range(n)] for i in range(n)]
    return eigenvalues_general(Acl)


def test_double_integrator():
    A = [[1.0, 1.0], [0.0, 1.0]]
    B = [[0.0], [1.0]]
    K = ackermann(A, B, [0.3, 0.5])
    eigs = sorted(e.real for e in _closed_loop_eigs(A, B, K))
    for got, want in zip(eigs, [0.3, 0.5]):
        assert close(got, want)


def test_scalar():
    K = ackermann([[2.0]], [[1.0]], [0.4])
    assert close(K[0][0], 1.6)


def test_three_state():
    A = [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-0.5, 0.3, 0.8]]
    B = [[0.0], [0.0], [1.0]]
    poles = [0.1, 0.2, 0.3]
    K = ackermann(A, B, poles)
    eigs = sorted(e.real for e in _closed_loop_eigs(A, B, K))
    for got, want in zip(eigs, poles):
        assert close(got, want, 1e-5)


def test_stabilizes_unstable_system():
    A = [[1.5, 0.5], [0.0, 1.2]]
    B = [[1.0], [1.0]]
    K = ackermann(A, B, [0.2, 0.4])
    for e in _closed_loop_eigs(A, B, K):
        assert abs(complex(e)) < 1.0


def test_errors():
    with pytest.raises(ValueError):
        ackermann([[1, 0], [0, 1]], [[1], [0]], [0.5])            # wrong pole count
    with pytest.raises(ValueError):
        ackermann([[1, 0], [0, 1]], [[1, 0], [0, 1]], [0.5, 0.5])  # multi-input
