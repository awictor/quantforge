"""Quadratic residues and multiplicative structure modulo a prime."""

import random

import pytest

from quantforge import (
    legendre_symbol,
    jacobi_symbol,
    tonelli_shanks,
    multiplicative_order,
    primitive_root,
)
from quantforge import mod_pow


def test_legendre_vs_brute():
    for p in (3, 5, 7, 11, 13, 17, 101):
        residues = {(x * x) % p for x in range(1, p)}
        for a in range(1, p):
            assert legendre_symbol(a, p) == (1 if a in residues else -1)


def test_jacobi():
    for p in (3, 5, 7, 11, 13, 17):
        for a in range(p):
            assert jacobi_symbol(a, p) == legendre_symbol(a, p)
    assert jacobi_symbol(2, 15) == jacobi_symbol(2, 3) * jacobi_symbol(2, 5)


def test_tonelli_shanks():
    rng = random.Random(1)
    for p in (7, 13, 17, 97, 101, 1009):
        for _ in range(50):
            x = rng.randint(1, p - 1)
            a = (x * x) % p
            r = tonelli_shanks(a, p)
            assert (r * r) % p == a
    p = 7
    nr = next(a for a in range(1, p) if legendre_symbol(a, p) == -1)
    with pytest.raises(ValueError):
        tonelli_shanks(nr, p)


def test_multiplicative_order():
    for p in (7, 11, 13, 17, 101):
        for a in range(2, p):
            o = multiplicative_order(a, p)
            assert (p - 1) % o == 0
            assert mod_pow(a, o, p) == 1


def test_primitive_root_generates_group():
    for p in (3, 5, 7, 11, 13, 17, 23, 101):
        g = primitive_root(p)
        assert multiplicative_order(g, p) == p - 1
        assert {mod_pow(g, k, p) for k in range(p - 1)} == set(range(1, p))


def test_validation():
    with pytest.raises(ValueError):
        legendre_symbol(3, 4)
    with pytest.raises(ValueError):
        jacobi_symbol(3, 4)
    with pytest.raises(ValueError):
        primitive_root(4)
