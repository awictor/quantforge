"""FX cross rates and triangular arbitrage."""

import pytest

from quantforge import cross_rate, triangular_arbitrage, is_arbitrage_free


def test_cross_rate_common_quote():
    # EURUSD 1.10, GBPUSD 1.25 -> EURGBP = 1.10 / 1.25.
    assert cross_rate(1.10, 1.25) == pytest.approx(0.88)


def test_cross_rate_chained():
    assert cross_rate(1.10, 0.9, via_is_quote=False) == pytest.approx(0.99)


def test_arbitrage_free_loop():
    eurusd, usdjpy = 1.10, 150.0
    jpyeur = 1 / (eurusd * usdjpy)
    assert triangular_arbitrage(eurusd, usdjpy, jpyeur) == pytest.approx(1.0)
    assert is_arbitrage_free(eurusd, usdjpy, jpyeur)


def test_mispriced_loop_signals_profit():
    assert triangular_arbitrage(1.10, 150, 1 / (1.05 * 150)) > 1
    assert not is_arbitrage_free(1.10, 150, 1 / (1.05 * 150))


def test_reverse_direction_below_one():
    assert triangular_arbitrage(1.10, 150, 1 / (1.15 * 150)) < 1


def test_validation():
    with pytest.raises(ValueError):
        cross_rate(0, 1.25)
    with pytest.raises(ValueError):
        triangular_arbitrage(1.1, -1, 0.5)
