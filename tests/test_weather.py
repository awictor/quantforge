"""Weather derivatives: degree days and temperature-index options."""

import math

import pytest

from quantforge import (
    heating_degree_days, cooling_degree_days, degree_day_index,
    degree_day_swap_payoff, degree_day_option,
    degree_day_swap_rate, degree_day_collar, degree_day_option_mc,
    seasonal_mean_temperature, expected_temperature, temperature_variance,
)


TEMPS = [60, 62, 70, 72, 55, 68, 64]


def test_hdd_cdd_complementary():
    # HDD - CDD = sum(base - T) exactly.
    diff = heating_degree_days(TEMPS, 65) - cooling_degree_days(TEMPS, 65)
    assert diff == pytest.approx(sum(65 - t for t in TEMPS))


def test_hdd_value():
    assert heating_degree_days(TEMPS, 65) == 5 + 3 + 10 + 1


def test_cdd_value():
    assert cooling_degree_days(TEMPS, 65) == 5 + 7 + 3


def test_degree_day_index_dispatch():
    assert degree_day_index(TEMPS, 65, "HDD") == heating_degree_days(TEMPS, 65)
    assert degree_day_index(TEMPS, 65, "CDD") == cooling_degree_days(TEMPS, 65)


def test_swap_payoff_linear():
    assert degree_day_swap_payoff(100, 80, 20, 1) == 400
    assert degree_day_swap_payoff(100, 80, 20, -1) == -400


def test_option_put_call_parity():
    c = degree_day_option(100, 90, 15, 0.03, 0.5, 20, True)
    p = degree_day_option(100, 90, 15, 0.03, 0.5, 20, False)
    assert c - p == pytest.approx(math.exp(-0.03 * 0.5) * 20 * (100 - 90), abs=1e-6)


def test_option_zero_vol_intrinsic():
    assert degree_day_option(100, 90, 0.0, 0.03, 0.5, 20, True) == pytest.approx(
        math.exp(-0.03 * 0.5) * 20 * 10, abs=1e-6)


def test_cap_reduces_call_value():
    uncapped = degree_day_option(100, 90, 15, 0.03, 0.5, 20, True)
    capped = degree_day_option(100, 90, 15, 0.03, 0.5, 20, True, cap=15)
    assert capped < uncapped
    assert capped <= math.exp(-0.03 * 0.5) * 20 * 15 + 1e-6


def test_seasonal_mean_sinusoid():
    assert seasonal_mean_temperature(0, 50, 0.001, 20, 80, 365) == pytest.approx(
        50 + 20 * math.sin(2 * math.pi * (-80) / 365))


def test_expected_temperature_zero_horizon():
    assert expected_temperature(70, 55, 55, 0.2, 0) == pytest.approx(70)


def test_expected_temperature_long_run():
    assert expected_temperature(70, 55, 60, 0.2, 200) == pytest.approx(60, abs=1e-6)


def test_expected_temperature_reverts():
    d1 = abs(expected_temperature(70, 55, 55, 0.2, 1) - 55)
    d2 = abs(expected_temperature(70, 55, 55, 0.2, 5) - 55)
    assert d2 < d1


def test_temperature_variance_limits():
    assert temperature_variance(3.0, 0.2, 0) == pytest.approx(0.0, abs=1e-12)
    assert temperature_variance(3.0, 0.2, 500) == pytest.approx(9 / (2 * 0.2), abs=1e-6)


def test_temperature_variance_monotone():
    vs = [temperature_variance(3.0, 0.2, h) for h in (1, 5, 10, 20, 50)]
    assert all(vs[i] < vs[i + 1] for i in range(len(vs) - 1))


def test_temperature_model_validation():
    with pytest.raises(ValueError):
        temperature_variance(3.0, 0, 1)
    with pytest.raises(ValueError):
        expected_temperature(70, 55, 55, -1, 1)


def test_swap_rate_zeroes_swap():
    ei = 120.0
    assert degree_day_swap_rate(ei) == ei
    assert degree_day_swap_payoff(ei, degree_day_swap_rate(ei), 20) == pytest.approx(0.0)


def test_collar_same_strike_is_forward():
    c = degree_day_collar(120, 100, 100, 15, 0.03, 0.5, 20)
    assert c == pytest.approx(math.exp(-0.03 * 0.5) * 20 * (120 - 100), abs=1e-6)


def test_collar_validation():
    with pytest.raises(ValueError):
        degree_day_collar(120, 100, 110, 15, 0.03, 0.5, 20)


def test_option_mc_deterministic():
    means = [75.0] * 90
    a = degree_day_option_mc(means, 3.0, 65.0, 880, 0.03, 0.5, 20, "CDD", True, 5000, 7)
    b = degree_day_option_mc(means, 3.0, 65.0, 880, 0.03, 0.5, 20, "CDD", True, 5000, 7)
    assert a == b


@pytest.mark.slow
def test_mc_validates_bachelier_option():
    means = [75.0] * 90
    ds, base = 3.0, 65.0
    ei = cooling_degree_days(means, base)
    sig_idx = ds * math.sqrt(90)
    bach = degree_day_option(ei, 880, sig_idx, 0.03, 0.5, 20, True)
    mcv = degree_day_option_mc(means, ds, base, 880, 0.03, 0.5, 20, "CDD", True,
                               80000, 7)
    assert abs(bach - mcv) / bach < 0.1


def test_validation():
    with pytest.raises(ValueError):
        degree_day_index(TEMPS, 65, "XYZ")
    with pytest.raises(ValueError):
        degree_day_option(100, 90, 15, 0.03, 0.5, 20, True, cap=-1)
    with pytest.raises(ValueError):
        degree_day_option(100, 90, -1, 0.03, 0.5, 20)
