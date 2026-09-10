"""Tests for whole-chain implied-vol smile inversion."""

import math

import pytest

from quantforge import implied_vol_smile, call_price, put_price, OptionType


def test_recovers_known_smile():
    S, t, r = 100.0, 1.0, 0.05
    strikes = [80, 90, 100, 110, 120]
    vols = [0.28, 0.24, 0.22, 0.225, 0.24]
    prices = [call_price(S, K, t, r, v) for K, v in zip(strikes, vols)]
    smile = implied_vol_smile(strikes, prices, S, t, r, OptionType.CALL)
    assert len(smile) == 5
    for (k, iv), v in zip(smile, vols):
        assert iv == pytest.approx(v, abs=1e-5)


def test_log_moneyness_axis():
    S, t, r = 100.0, 1.0, 0.0   # forward = S
    strikes = [90, 100, 110]
    prices = [call_price(S, K, t, r, 0.2) for K in strikes]
    smile = implied_vol_smile(strikes, prices, S, t, r)
    ks = [k for k, _ in smile]
    assert ks[0] == pytest.approx(math.log(90 / 100))
    assert ks[1] == pytest.approx(0.0)
    assert ks[2] == pytest.approx(math.log(110 / 100))


def test_drops_arbitrage_quotes():
    S, t, r = 100.0, 1.0, 0.05
    # 200 exceeds the forward (arb) so that quote is dropped; 5 is fine.
    smile = implied_vol_smile([100, 110], [200, 5], S, t, r, OptionType.CALL)
    assert len(smile) == 1


def test_sorted_by_strike():
    S, t, r = 100.0, 1.0, 0.05
    strikes = [110, 90, 100]   # unsorted
    prices = [call_price(S, K, t, r, 0.2) for K in strikes]
    smile = implied_vol_smile(strikes, prices, S, t, r)
    ks = [k for k, _ in smile]
    assert ks == sorted(ks)


def test_put_smile():
    S, t, r = 100.0, 0.5, 0.03
    strikes = [90, 100, 110]
    vols = [0.26, 0.22, 0.2]
    prices = [put_price(S, K, t, r, v) for K, v in zip(strikes, vols)]
    smile = implied_vol_smile(strikes, prices, S, t, r, OptionType.PUT)
    for (_, iv), v in zip(smile, vols):
        assert iv == pytest.approx(v, abs=1e-5)


def test_length_mismatch_rejected():
    with pytest.raises(ValueError):
        implied_vol_smile([90, 100], [1.0], 100, 1.0, 0.05)
