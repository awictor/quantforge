"""Tests for Merton jump-diffusion pricing."""

import math

import pytest

from quantforge import merton_jump_price as mj, call_price, put_price, OptionType


@pytest.mark.parametrize("K", [80, 100, 120])
@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_zero_intensity_reduces_to_bsm(K, ot):
    # lam = 0: no jumps, so the series is a single Black-Scholes term.
    v = mj(100, K, 1.0, 0.05, 0.2, lam=0.0, mu_j=0.0, sigma_j=0.0, option_type=ot)
    ref = (call_price(100, K, 1.0, 0.05, 0.2) if ot is OptionType.CALL
           else put_price(100, K, 1.0, 0.05, 0.2))
    assert v == pytest.approx(ref, abs=1e-9)


def test_put_call_parity():
    kw = dict(S=100, K=95, t=1.0, r=0.05, sigma=0.2, lam=1.0, mu_j=-0.1, sigma_j=0.15)
    c = mj(**kw, option_type=OptionType.CALL)
    p = mj(**kw, option_type=OptionType.PUT)
    # Correct martingale forward => parity holds exactly.
    assert c - p == pytest.approx(100 - 95 * math.exp(-0.05), abs=1e-8)


def test_reference_value_matches_monte_carlo():
    # Cross-checked against a compensated-drift Monte Carlo (CF 12.761, MC ~12.75).
    v = mj(100, 100, 1.0, 0.05, 0.2, lam=1.0, mu_j=-0.1, sigma_j=0.15,
           option_type=OptionType.CALL)
    assert v == pytest.approx(12.761, abs=0.05)


def test_jumps_add_value_to_atm_option():
    # Adding jump risk raises the value of an ATM option (fatter tails).
    base = call_price(100, 100, 1.0, 0.05, 0.2)
    with_jumps = mj(100, 100, 1.0, 0.05, 0.2, lam=2.0, mu_j=0.0, sigma_j=0.2)
    assert with_jumps > base


def test_series_converges_for_high_intensity():
    # Large lam*t must still sum to a finite, positive price.
    v = mj(100, 100, 1.0, 0.05, 0.2, lam=20.0, mu_j=-0.05, sigma_j=0.1,
           option_type=OptionType.CALL)
    assert math.isfinite(v) and v > 0


def test_negative_jump_mean_skews_puts_up():
    # Negative mean jump (crash risk) makes OTM puts richer than the pure
    # diffusion value.
    otm_put_bsm = put_price(100, 85, 1.0, 0.05, 0.2)
    otm_put_jump = mj(100, 85, 1.0, 0.05, 0.2, lam=1.0, mu_j=-0.3, sigma_j=0.2,
                      option_type=OptionType.PUT)
    assert otm_put_jump > otm_put_bsm


def test_invalid_inputs_rejected():
    with pytest.raises(ValueError):
        mj(100, 100, 1.0, 0.05, 0.2, lam=-1.0, mu_j=0.0, sigma_j=0.1)
    with pytest.raises(ValueError):
        mj(100, 100, 1.0, 0.05, 0.2, lam=1.0, mu_j=0.0, sigma_j=-0.1)


def test_zero_time_is_intrinsic():
    assert mj(110, 100, 0.0, 0.05, 0.2, lam=1.0, mu_j=0.0, sigma_j=0.1,
              option_type=OptionType.CALL) == pytest.approx(10.0)
