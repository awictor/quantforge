"""Pade approximants and continued-fraction evaluation."""

import math

import pytest

from quantforge import pade, pade_eval, lentz_continued_fraction


def test_exp_pade_2_2_coefficients():
    coeffs = [1 / math.factorial(k) for k in range(5)]
    num, den = pade(coeffs, 2, 2)
    assert abs(num[0] - 1.0) < 1e-12 and abs(num[1] - 0.5) < 1e-12
    assert abs(num[2] - 1 / 12) < 1e-12
    assert abs(den[1] + 0.5) < 1e-12 and abs(den[2] - 1 / 12) < 1e-12


def test_pade_beats_taylor():
    coeffs = [1 / math.factorial(k) for k in range(5)]
    num, den = pade(coeffs, 2, 2)
    p = pade_eval(num, den, 1.0)
    t = sum(coeffs[k] for k in range(5))
    assert abs(p - math.e) < abs(t - math.e)


def test_pade_series_matches_through_order():
    coeffs = [1 / math.factorial(k) for k in range(5)]
    num, den = pade(coeffs, 2, 2)

    def series(num, den, order):
        c = []
        for k in range(order + 1):
            ck = (num[k] if k < len(num) else 0.0)
            for j in range(1, min(k, len(den) - 1) + 1):
                ck -= den[j] * c[k - j]
            c.append(ck / den[0])
        return c

    ps = series(num, den, 4)
    assert all(abs(ps[k] - coeffs[k]) < 1e-9 for k in range(5))


def test_geometric_exact():
    g = [(-1) ** k for k in range(5)]        # 1/(1+x)
    num, den = pade(g, 0, 1)
    assert abs(pade_eval(num, den, 3.0) - 0.25) < 1e-12


def test_log_pade_converges_beyond_taylor_radius():
    lc = [0.0] + [(-1) ** (k + 1) / k for k in range(1, 9)]   # ln(1+x)
    num, den = pade(lc, 4, 4)
    assert abs(pade_eval(num, den, 2.0) - math.log(3)) < 1e-3   # Taylor diverges at x=2


def test_lentz_golden_ratio():
    phi = lentz_continued_fraction(lambda k: 1.0, lambda k: 1.0)
    assert abs(phi - (1 + math.sqrt(5)) / 2) < 1e-12


def test_lentz_tangent():
    x = 1.0
    a = lambda k: x if k == 1 else -x * x
    b = lambda k: 0.0 if k == 0 else (2 * k - 1)
    assert abs(lentz_continued_fraction(a, b) - math.tan(1.0)) < 1e-10


def test_validation():
    with pytest.raises(ValueError):
        pade([1.0, 1.0], 2, 2)                    # not enough coeffs
    with pytest.raises(ValueError):
        pade([1.0, 1.0, 1.0], -1, 1)
