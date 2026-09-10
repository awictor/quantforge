"""Tests for quanto options."""

import pytest

from quantforge import quanto_option, price, OptionType


BASE = dict(S=100, K=100, t=1.0, r_domestic=0.03, r_foreign=0.05,
            sigma_asset=0.25, sigma_fx=0.1)


def test_zero_correlation_is_foreign_carry_bsm():
    # rho = 0: no quanto adjustment; carry = r_foreign, discount = r_domestic.
    q = quanto_option(**BASE, rho=0.0, option_type=OptionType.CALL)
    ref = price(100, 100, 1.0, 0.03, 0.25, OptionType.CALL, b=0.05)
    assert q == pytest.approx(ref, abs=1e-9)


def test_positive_correlation_lowers_call():
    base = quanto_option(**BASE, rho=0.0, option_type=OptionType.CALL)
    pos = quanto_option(**BASE, rho=0.5, option_type=OptionType.CALL)
    assert pos < base


def test_negative_correlation_raises_call():
    base = quanto_option(**BASE, rho=0.0, option_type=OptionType.CALL)
    neg = quanto_option(**BASE, rho=-0.5, option_type=OptionType.CALL)
    assert neg > base


def test_matches_monte_carlo():
    # Cross-checked against a quanto-drift Monte Carlo (~11.80 at rho=0.5).
    q = quanto_option(**BASE, rho=0.5, option_type=OptionType.CALL)
    assert q == pytest.approx(11.80, abs=0.1)


def test_put_call_parity_with_quanto_carry():
    import math
    q_c = quanto_option(**BASE, rho=0.3, option_type=OptionType.CALL)
    q_p = quanto_option(**BASE, rho=0.3, option_type=OptionType.PUT)
    bq = 0.05 - 0.3 * 0.25 * 0.1
    lhs = q_c - q_p
    rhs = 100 * math.exp((bq - 0.03) * 1.0) - 100 * math.exp(-0.03 * 1.0)
    assert lhs == pytest.approx(rhs, abs=1e-9)


def test_zero_fx_vol_removes_adjustment():
    # sigma_fx = 0: correlation term vanishes regardless of rho.
    a = quanto_option(**{**BASE, "sigma_fx": 0.0}, rho=0.8, option_type=OptionType.CALL)
    b = quanto_option(**{**BASE, "sigma_fx": 0.0}, rho=0.0, option_type=OptionType.CALL)
    assert a == pytest.approx(b, abs=1e-12)


def test_rejects_bad_correlation():
    with pytest.raises(ValueError):
        quanto_option(**BASE, rho=1.5)
