"""Black-Litterman optimal weights (portopt module)."""

import pytest

from quantforge import black_litterman_weights


COV = [[0.04, 0.01, 0.0], [0.01, 0.09, 0.02], [0.0, 0.02, 0.16]]
MW = [0.5, 0.3, 0.2]


def test_no_views_recovers_market_weights():
    w = black_litterman_weights(COV, MW, [], [])
    assert w == pytest.approx(MW, abs=1e-9)


def test_no_views_raw_weights_are_market():
    raw = black_litterman_weights(COV, MW, [], [], normalize=False)
    assert raw == pytest.approx(MW, abs=1e-9)


def test_bullish_view_tilts_toward_asset():
    w = black_litterman_weights(COV, MW, [[1.0, 0.0, 0.0]], [0.15])
    assert w[0] > MW[0]
    assert sum(w) == pytest.approx(1.0, abs=1e-9)


def test_bearish_view_tilts_away():
    w = black_litterman_weights(COV, MW, [[1.0, 0.0, 0.0]], [-0.05])
    assert w[0] < MW[0]


def test_normalized_sums_to_one():
    w = black_litterman_weights(COV, MW, [[0.0, 1.0, -1.0]], [0.03])
    assert sum(w) == pytest.approx(1.0, abs=1e-9)
