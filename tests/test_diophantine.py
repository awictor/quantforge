"""Diophantine equations: linear solutions, sqrt continued fractions, Pell."""

import math
import random

import pytest

from quantforge import linear_diophantine, sqrt_continued_fraction, pell_fundamental


def test_linear_diophantine():
    rng = random.Random(1)
    for _ in range(2000):
        a = rng.randint(-50, 50)
        b = rng.randint(-50, 50)
        if a == 0 and b == 0:
            continue
        c = rng.randint(-100, 100)
        res = linear_diophantine(a, b, c)
        if c % math.gcd(a, b) != 0:
            assert res is None
            continue
        x0, y0, dx, dy = res
        assert a * x0 + b * y0 == c
        assert a * (x0 + dx) + b * (y0 + dy) == c
    assert linear_diophantine(2, 4, 5) is None


def test_sqrt_continued_fraction():
    assert sqrt_continued_fraction(2) == (1, [2])
    assert sqrt_continued_fraction(7) == (2, [1, 1, 1, 4])
    assert sqrt_continued_fraction(23) == (4, [1, 3, 1, 8])
    for n in (2, 3, 5, 7, 13, 23, 61):
        a0, per = sqrt_continued_fraction(n)
        terms = [a0] + per * 20
        h_prev, h = 1, terms[0]
        k_prev, k = 0, 1
        for a in terms[1:]:
            h_prev, h = h, a * h + h_prev
            k_prev, k = k, a * k + k_prev
        assert abs(h / k - math.sqrt(n)) < 1e-9


def test_pell_fundamental():
    known = {2: (3, 2), 3: (2, 1), 5: (9, 4), 6: (5, 2), 7: (8, 3),
             13: (649, 180), 61: (1766319049, 226153980)}
    for n, sol in known.items():
        assert pell_fundamental(n) == sol
        x, y = sol
        assert x * x - n * y * y == 1


def test_validation():
    with pytest.raises(ValueError):
        sqrt_continued_fraction(9)                    # perfect square
    with pytest.raises(ValueError):
        linear_diophantine(0, 0, 5)
    with pytest.raises(ValueError):
        pell_fundamental(16)                          # perfect square
