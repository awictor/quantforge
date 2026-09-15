import random

from quantforge import lanczos_eigenvalues
from quantforge.pca import jacobi_eigen


def close(a, b, tol=1e-4):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _rand_sym(n, seed=0):
    random.seed(seed)
    M = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            M[i][j] = M[j][i] = random.gauss(0, 1)
    return M


def test_diagonal_matrix():
    diag = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0]
    D = [[diag[i] if i == j else 0.0 for j in range(6)] for i in range(6)]
    rv = sorted(lanczos_eigenvalues(D, 6, 6))
    for a, b in zip(rv, sorted(diag)):
        assert close(a, b, 1e-6)


def test_extremal_ritz_values():
    M = _rand_sym(30)
    vals = sorted(jacobi_eigen(M)[0])
    rv = sorted(lanczos_eigenvalues(M, 30, 20))
    assert close(rv[0], vals[0], 1e-3)
    assert close(rv[-1], vals[-1], 1e-3)


def test_full_spectrum_when_m_equals_n():
    M = _rand_sym(30)
    vals = sorted(jacobi_eigen(M)[0])
    rv = sorted(lanczos_eigenvalues(M, 30, 30))
    assert max(abs(rv[i] - vals[i]) for i in range(30)) < 1e-3


def test_matvec_callable():
    M = _rand_sym(30)
    vals = sorted(jacobi_eigen(M)[0])
    mv = lambda v: [sum(M[i][j] * v[j] for j in range(30)) for i in range(30)]
    rv = sorted(lanczos_eigenvalues(mv, 30, 20))
    assert close(rv[0], vals[0], 1e-3)
    assert close(rv[-1], vals[-1], 1e-3)


def test_spd_largest_eigenvalue():
    M = _rand_sym(20)
    n = 20
    A = [[sum(M[k][i] * M[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    av = sorted(jacobi_eigen(A)[0])
    rv = sorted(lanczos_eigenvalues(A, n, 15))
    assert close(rv[-1], av[-1], 1e-3)
    assert rv[0] > -1e-6
