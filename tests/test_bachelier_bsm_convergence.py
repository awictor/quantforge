"""Bachelier and Black-Scholes agree in the low-vol ATM limit.

At the money, matching the normal vol to the lognormal vol via
``sigma_normal = sigma_lognormal * F`` makes the two models coincide as the
volatility shrinks (both approach the same small-move Gaussian), with the gap
growing smoothly in vol.
"""

import pytest

from quantforge import bachelier_price, call_price, put_price, OptionType


def test_atm_convergence_shrinks_with_vol():
    S, K, t, r = 100.0, 100.0, 1.0, 0.0
    diffs = []
    for sig in (0.05, 0.1, 0.2):
        bs = call_price(S, K, t, r, sig)
        ba = bachelier_price(S, K, t, r, sig * S)   # sigma_normal = sigma_ln * F
        diffs.append(abs(bs - ba))
    # Monotonically increasing gap as vol grows; tiny at 5% vol.
    assert diffs[0] < 1e-3
    assert diffs[0] < diffs[1] < diffs[2]


def test_low_vol_call_and_put_match():
    S, K, t, r, sig = 100.0, 100.0, 1.0, 0.0, 0.05
    assert bachelier_price(S, K, t, r, sig * S, OptionType.CALL) == pytest.approx(
        call_price(S, K, t, r, sig), abs=1e-3)
    assert bachelier_price(S, K, t, r, sig * S, OptionType.PUT) == pytest.approx(
        put_price(S, K, t, r, sig), abs=1e-3)
