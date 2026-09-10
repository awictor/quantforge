"""Richardson-extrapolated (Broadie-Detemple) Leisen-Reimer American price."""

import pytest

from quantforge import (
    OptionType,
    leisen_reimer_american_accel as accel,
    leisen_reimer_price as lr,
    american_price as crr,
    call_price,
)


AMERICAN_CASES = [
    (42, 40, 0.75, 0.04, 0.35, OptionType.CALL, -0.04),
    (100, 100, 0.5, 0.05, 0.30, OptionType.PUT, 0.05),
    (110, 100, 0.5, 0.10, 0.25, OptionType.PUT, 0.10),
    (90, 100, 1.0, 0.08, 0.20, OptionType.PUT, 0.08),
]


@pytest.mark.slow
@pytest.mark.parametrize("S,K,t,r,sigma,ot,b", AMERICAN_CASES)
def test_accel_more_accurate_than_plain(S, K, t, r, sigma, ot, b):
    ref = crr(S, K, t, r, sigma, ot, b=b, steps=6000)
    plain = lr(S, K, t, r, sigma, ot, b=b, steps=101, american=True)
    ac = accel(S, K, t, r, sigma, ot, b=b, steps=101)
    # Extrapolation should not be worse, and reaches ~a few mils.
    assert abs(ac - ref) <= abs(plain - ref) + 1e-9
    assert abs(ac - ref) < 3e-3


@pytest.mark.slow
@pytest.mark.parametrize("S,K,t,r,sigma,ot,b", AMERICAN_CASES)
def test_accel_matches_converged_crr(S, K, t, r, sigma, ot, b):
    ref = crr(S, K, t, r, sigma, ot, b=b, steps=6000)
    assert accel(S, K, t, r, sigma, ot, b=b, steps=151) == pytest.approx(ref,
                                                                         abs=3e-3)


def test_no_dividend_call_equals_european():
    # b = r: an American call is never exercised early -> European value.
    ac = accel(100, 100, 1.0, 0.05, 0.2, OptionType.CALL)
    assert ac == pytest.approx(call_price(100, 100, 1.0, 0.05, 0.2), abs=2e-3)


def test_accel_above_european_put():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    am = accel(S, K, t, r, sigma, OptionType.PUT)
    eu = lr(S, K, t, r, sigma, OptionType.PUT, steps=201, american=False)
    assert am > eu


def test_deep_itm_put_dominates_intrinsic():
    S, K, t, r, sigma = 60, 100, 0.5, 0.05, 0.3
    assert accel(S, K, t, r, sigma, OptionType.PUT) >= (K - S) - 1e-6
