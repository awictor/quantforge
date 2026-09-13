"""Nonparametric survival estimators: Kaplan-Meier, Nelson-Aalen."""

import math

import pytest

from quantforge import kaplan_meier, nelson_aalen, survival_at


def test_km_no_censoring_is_one_minus_ecdf():
    t = [1, 2, 3, 4, 5]
    e = [1, 1, 1, 1, 1]
    _, s = kaplan_meier(t, e)
    assert [round(x, 6) for x in s] == [0.8, 0.6, 0.4, 0.2, 0.0]


def test_km_non_increasing_in_unit_interval():
    t = [1, 2, 2, 3, 5, 8]
    e = [1, 1, 0, 1, 1, 0]
    _, s = kaplan_meier(t, e)
    assert all(s[i] >= s[i + 1] for i in range(len(s) - 1))
    assert all(0.0 <= x <= 1.0 for x in s)


def test_km_censored_example_by_hand():
    t = [2, 3, 3, 5, 7]
    e = [1, 1, 0, 1, 1]      # the second time-3 observation is censored
    ts, s = kaplan_meier(t, e)
    assert ts == [2, 3, 5, 7]
    assert abs(s[0] - 0.8) < 1e-9      # 1 - 1/5
    assert abs(s[1] - 0.6) < 1e-9      # 0.8 * (1 - 1/4)


def test_nelson_aalen_cumulative():
    t = [2, 3, 3, 5, 7]
    e = [1, 1, 0, 1, 1]
    _, h = nelson_aalen(t, e)
    # H = 1/5, +1/4, +1/2, +1/1
    assert abs(h[0] - 0.2) < 1e-9
    assert abs(h[1] - 0.45) < 1e-9
    assert abs(h[-1] - 1.95) < 1e-9


def test_na_survival_close_to_km_early():
    t = [2, 3, 5, 8, 10, 12]
    e = [1, 1, 1, 1, 1, 1]
    _, s_km = kaplan_meier(t, e)
    _, h = nelson_aalen(t, e)
    # exp(-H) tracks KM when the risk set is large (first step).
    assert abs(math.exp(-h[0]) - s_km[0]) < 0.02


def test_survival_at_step():
    t = [2, 3, 3, 5, 7]
    e = [1, 1, 0, 1, 1]
    assert abs(survival_at(t, e, 4) - 0.6) < 1e-9    # no event in (3, 4]
    assert survival_at(t, e, 1) == 1.0               # before the first event


def test_validation():
    with pytest.raises(ValueError):
        kaplan_meier([], [])
    with pytest.raises(ValueError):
        nelson_aalen([1, 2], [1])
    with pytest.raises(ValueError):
        kaplan_meier([-1, 2], [1, 1])
    with pytest.raises(ValueError):
        survival_at([1, 2], [1, 1], 1.5, estimator="bad")
