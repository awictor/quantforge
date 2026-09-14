"""Association strength for categorical contingency tables."""

import math
import random

import pytest

from quantforge import (
    cramers_v,
    phi_coefficient,
    tschuprow_t,
    contingency_coefficient,
)


def test_independence_is_zero():
    assert cramers_v([[10, 20], [20, 40]]) < 1e-9      # proportional rows
    assert phi_coefficient([[10, 20], [20, 40]]) < 1e-9


def test_perfect_association():
    assert abs(cramers_v([[50, 0], [0, 50]]) - 1.0) < 1e-9


def test_phi_equals_cramers_v_and_closed_form():
    rng = random.Random(1)
    for _ in range(500):
        a, b, c, d = (rng.randint(1, 50) for _ in range(4))
        t = [[a, b], [c, d]]
        assert abs(phi_coefficient(t) - cramers_v(t)) < 1e-12
        phi = abs(a * d - b * c) / math.sqrt((a + b) * (c + d) * (a + c) * (b + d))
        assert abs(phi_coefficient(t) - phi) < 1e-9


def test_square_table_v_equals_tschuprow():
    t = [[20, 5, 5], [5, 20, 5], [5, 5, 20]]
    assert abs(cramers_v(t) - 0.5) < 1e-9
    assert abs(cramers_v(t) - tschuprow_t(t)) < 1e-9    # equal for square tables
    assert 0.0 < contingency_coefficient(t) < 1.0


def test_validation():
    with pytest.raises(ValueError):
        phi_coefficient([[1, 2, 3], [4, 5, 6]])         # not 2x2
    with pytest.raises(ValueError):
        cramers_v([[0, 0], [0, 0]])                     # empty
