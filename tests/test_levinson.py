"""Tests for Levinson-Durbin: Toeplitz solves and AR fitting, vs dense linear algebra."""

import random

import pytest

from quantforge.levinson import solve_toeplitz, levinson_durbin


def _dense_solve(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c]:
                f = M[r][c] / M[c][c]
                for k in range(c, n + 1):
                    M[r][k] -= f * M[c][k]
    return [M[i][n] / M[i][i] for i in range(n)]


def _close(a, b, tol=1e-6):
    return len(a) == len(b) and all(abs(x - y) <= tol * max(1, abs(x), abs(y)) for x, y in zip(a, b))


def test_fuzz_toeplitz_vs_dense():
    rng = random.Random(461)
    for _ in range(3000):
        n = rng.randint(1, 8)
        r = [rng.uniform(5, 10)] + [rng.uniform(-2, 2) for _ in range(n - 1)]
        b = [rng.uniform(-5, 5) for _ in range(n)]
        T = [[r[abs(i - j)] for j in range(n)] for i in range(n)]
        assert _close(solve_toeplitz(r, b), _dense_solve(T, b))


def test_ar_coeffs_match_yule_walker():
    r = [1.0, 0.5, 0.3, 0.1]
    ar, err, refl = levinson_durbin(r)
    p = len(r) - 1
    Rt = [[r[abs(i - j)] for j in range(p)] for i in range(p)]
    rhs = [r[i + 1] for i in range(p)]
    assert _close(ar, _dense_solve(Rt, rhs))
    assert err > 0
    assert all(abs(k) < 1 for k in refl)


def test_ar1_exact_recovery():
    rho = 0.7
    r = [rho ** k for k in range(5)]
    ar, err, refl = levinson_durbin(r)
    assert _close([ar[0]], [rho], 1e-6)
    assert _close(ar[1:], [0, 0, 0], 1e-6)


def test_toeplitz_identity():
    # diagonal Toeplitz -> x = b / r0
    assert _close(solve_toeplitz([2, 0, 0], [4, 6, 8]), [2, 3, 4])


def test_solution_satisfies_system():
    rng = random.Random(462)
    n = 6
    r = [rng.uniform(4, 8)] + [rng.uniform(-1, 1) for _ in range(n - 1)]
    b = [rng.uniform(-3, 3) for _ in range(n)]
    x = solve_toeplitz(r, b)
    for i in range(n):
        val = sum(r[abs(i - j)] * x[j] for j in range(n))
        assert abs(val - b[i]) < 1e-8


def test_single_element():
    assert _close(solve_toeplitz([5], [10]), [2.0])


def test_empty():
    assert solve_toeplitz([], []) == []


def test_reflection_coeff_count():
    _, _, refl = levinson_durbin([1.0, 0.4, 0.2, 0.1, 0.05])
    assert len(refl) == 4


def test_zero_diagonal_raises():
    with pytest.raises(ValueError):
        solve_toeplitz([0, 1], [1, 1])
    with pytest.raises(ValueError):
        levinson_durbin([0, 1])


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        solve_toeplitz([1, 2], [1])


def test_reflection_white_noise():
    # white noise: r = [1, 0, 0, ...] -> all reflection coeffs 0, error 1
    ar, err, refl = levinson_durbin([1.0, 0.0, 0.0])
    assert _close(ar, [0.0, 0.0], 1e-9)
    assert abs(err - 1.0) < 1e-9
    assert all(abs(k) < 1e-12 for k in refl)
