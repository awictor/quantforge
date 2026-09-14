"""General (non-symmetric) eigenvalues via Faddeev-LeVerrier."""

import random

import pytest

from quantforge import (
    characteristic_polynomial,
    eigenvalues_general,
    determinant_from_charpoly,
)
from quantforge import jacobi_eigen, determinant


def test_characteristic_polynomial_known():
    assert characteristic_polynomial([[2, 0], [0, 3]]) == [1.0, -5.0, 6.0]


def test_eigenvalues_diagonal():
    evs = sorted(eigenvalues_general([[2, 0], [0, 3]]))
    assert abs(evs[0] - 2) < 1e-9 and abs(evs[1] - 3) < 1e-9


def test_complex_conjugate_pair():
    evs = eigenvalues_general([[0, -1], [1, 0]])         # rotation -> ±i
    assert all(isinstance(e, complex) for e in evs)
    assert abs(evs[0] - (-1j)) < 1e-9 and abs(evs[1] - 1j) < 1e-9


def test_matches_jacobi_on_symmetric():
    rng = random.Random(1)
    for _ in range(100):
        n = rng.randint(2, 4)
        M = [[rng.uniform(-3, 3) for _ in range(n)] for _ in range(n)]
        S = [[M[i][j] + M[j][i] for j in range(n)] for i in range(n)]
        eg = sorted(eigenvalues_general(S))
        jv = sorted(jacobi_eigen(S)[0])
        for a, b in zip(eg, jv):
            assert abs(a - b) < 1e-6


def test_sum_is_trace_product_is_det():
    rng = random.Random(2)
    for _ in range(100):
        n = rng.randint(2, 4)
        A = [[rng.uniform(-3, 3) for _ in range(n)] for _ in range(n)]
        evs = eigenvalues_general(A)
        s = sum(complex(e) for e in evs)
        assert abs(s.real - sum(A[i][i] for i in range(n))) < 1e-6 and abs(s.imag) < 1e-6
        prod = 1
        for e in evs:
            prod *= complex(e)
        assert abs(prod.real - determinant(A)) < 1e-5 and abs(prod.imag) < 1e-5


def test_determinant_from_charpoly():
    rng = random.Random(3)
    for _ in range(50):
        n = rng.randint(2, 4)
        A = [[rng.uniform(-3, 3) for _ in range(n)] for _ in range(n)]
        assert abs(determinant_from_charpoly(A) - determinant(A)) < 1e-6


def test_validation():
    with pytest.raises(ValueError):
        characteristic_polynomial([[1, 2]])              # non-square
