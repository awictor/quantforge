"""Property test: barrier in-out parity across all kinds/params.

For a zero-rebate single-barrier option, knock-in + knock-out must equal the
vanilla option at the same strike, for both calls and puts and every barrier
direction, regardless of where the barrier sits relative to spot and strike.
This is an exact analytic identity (the two events partition every path), so it
holds to machine precision for the Reiner-Rubinstein closed forms.
"""

import itertools

import pytest

from quantforge import barrier_option, price, Barrier, OptionType


IN_OUT = [(Barrier.DOWN_IN, Barrier.DOWN_OUT), (Barrier.UP_IN, Barrier.UP_OUT)]
STRIKES = [90, 100, 110]
BARRIERS = [80, 95, 105, 120]
OTS = [OptionType.CALL, OptionType.PUT]


@pytest.mark.parametrize("ki,ko", IN_OUT)
@pytest.mark.parametrize("K", STRIKES)
@pytest.mark.parametrize("H", BARRIERS)
@pytest.mark.parametrize("ot", OTS)
def test_in_out_parity(ki, ko, K, H, ot):
    S, t, r, sigma, b = 100.0, 1.0, 0.05, 0.25, 0.05
    v_in = barrier_option(S, K, H, t, r, sigma, ot, ki, b=b)
    v_out = barrier_option(S, K, H, t, r, sigma, ot, ko, b=b)
    vanilla = price(S, K, t, r, sigma, ot, b=b)
    assert v_in + v_out == pytest.approx(vanilla, abs=1e-8)


@pytest.mark.parametrize("ki,ko", IN_OUT)
def test_parity_holds_with_dividends(ki, ko):
    S, K, H, t, r, sigma, b = 100.0, 100.0, 90.0, 1.0, 0.06, 0.3, 0.06 - 0.03
    for ot in OTS:
        v_in = barrier_option(S, K, H, t, r, sigma, ot, ki, b=b)
        v_out = barrier_option(S, K, H, t, r, sigma, ot, ko, b=b)
        vanilla = price(S, K, t, r, sigma, ot, b=b)
        assert v_in + v_out == pytest.approx(vanilla, abs=1e-8)


def test_knock_out_never_exceeds_vanilla():
    # With no rebate a knock-out is worth at most the vanilla (parity + KI >= 0).
    S, t, r, sigma = 100.0, 1.0, 0.05, 0.25
    for K, H, ot, (_, ko) in itertools.product(STRIKES, BARRIERS, OTS, IN_OUT):
        v_out = barrier_option(S, K, H, t, r, sigma, ot, ko)
        vanilla = price(S, K, t, r, sigma, ot)
        assert v_out <= vanilla + 1e-8
