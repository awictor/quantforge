import random

from quantforge import (
    matrix_bandwidth,
    reverse_cuthill_mckee,
    apply_permutation,
    permute_vector,
    inverse_permutation,
)
from quantforge.lu import lu_solve


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _shuffled_path(n, seed=0):
    random.seed(seed)
    labels = list(range(n))
    random.shuffle(labels)
    A = [[0.0] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = 2.0
    for i in range(n - 1):
        a, b = labels[i], labels[i + 1]
        A[a][b] = A[b][a] = -1.0
    return A


def test_reduces_bandwidth_of_shuffled_path():
    A = _shuffled_path(20)
    bw0 = matrix_bandwidth(A)
    B = apply_permutation(A, reverse_cuthill_mckee(A))
    bw1 = matrix_bandwidth(B)
    assert bw1 < bw0
    assert bw1 <= 2


def test_permutation_is_valid():
    A = _shuffled_path(20)
    perm = reverse_cuthill_mckee(A)
    assert sorted(perm) == list(range(20))


def test_solve_consistency():
    A = _shuffled_path(20)
    n = 20
    for i in range(n):
        A[i][i] = 5.0
    b = [random.gauss(0, 1) for _ in range(n)]
    xref = lu_solve(A, b)
    perm = reverse_cuthill_mckee(A)
    xp = lu_solve(apply_permutation(A, perm), permute_vector(b, perm))
    inv = inverse_permutation(perm)
    xback = [xp[inv[i]] for i in range(n)]
    for a, c in zip(xback, xref):
        assert close(a, c)


def test_inverse_permutation():
    A = _shuffled_path(20)
    perm = reverse_cuthill_mckee(A)
    inv = inverse_permutation(perm)
    assert all(perm[inv[i]] == i for i in range(20))


def test_grid_bandwidth_not_worse():
    m = 5
    N = m * m
    G = [[0.0] * N for _ in range(N)]

    def idx(r, c):
        return r * m + c

    for r in range(m):
        for c in range(m):
            G[idx(r, c)][idx(r, c)] = 4.0
            if r + 1 < m:
                G[idx(r, c)][idx(r + 1, c)] = G[idx(r + 1, c)][idx(r, c)] = -1.0
            if c + 1 < m:
                G[idx(r, c)][idx(r, c + 1)] = G[idx(r, c + 1)][idx(r, c)] = -1.0
    bw0 = matrix_bandwidth(G)
    bw1 = matrix_bandwidth(apply_permutation(G, reverse_cuthill_mckee(G)))
    assert bw1 <= bw0
