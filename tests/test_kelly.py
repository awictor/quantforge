"""Tests for Kelly-criterion position sizing."""

import pytest

from quantforge import (
    kelly_fraction_binary, kelly_fraction_continuous, kelly_growth_rate,
)


def test_even_money_kelly_is_two_p_minus_one():
    # Even-money bet (b=1): f = 2p - 1.
    assert kelly_fraction_binary(0.6, 1.0) == pytest.approx(0.2)
    assert kelly_fraction_binary(0.75, 1.0) == pytest.approx(0.5)


def test_no_edge_is_zero():
    assert kelly_fraction_binary(0.5, 1.0) == pytest.approx(0.0)


def test_negative_edge_gives_negative_fraction():
    # A losing bet -> negative Kelly (caller clamps to zero / no bet).
    assert kelly_fraction_binary(0.4, 1.0) < 0


def test_favorable_odds_raise_fraction():
    # Higher payoff odds at the same win prob -> larger stake.
    lo = kelly_fraction_binary(0.5, 1.0)
    hi = kelly_fraction_binary(0.5, 3.0)
    assert hi > lo


def test_continuous_kelly_is_mu_over_variance():
    assert kelly_fraction_continuous(0.08, 0.04) == pytest.approx(2.0)


def test_fractional_kelly_scales():
    full = kelly_fraction_continuous(0.08, 0.04)
    half = kelly_fraction_continuous(0.08, 0.04, fraction=0.5)
    assert half == pytest.approx(0.5 * full)


def test_growth_maximized_at_full_kelly():
    mu, var = 0.08, 0.04
    f_star = kelly_fraction_continuous(mu, var)
    g_star = kelly_growth_rate(mu, var, f_star)
    for f in (0.7 * f_star, 0.9 * f_star, 1.1 * f_star, 1.3 * f_star):
        assert kelly_growth_rate(mu, var, f) <= g_star + 1e-12


def test_rejects_bad_inputs():
    with pytest.raises(ValueError):
        kelly_fraction_binary(1.5, 1.0)
    with pytest.raises(ValueError):
        kelly_fraction_binary(0.6, -1.0)
    with pytest.raises(ValueError):
        kelly_fraction_continuous(0.05, 0.0)
