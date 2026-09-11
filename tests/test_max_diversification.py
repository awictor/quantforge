"""Most-diversified portfolio (portopt module)."""

import pytest

from quantforge import (
    max_diversification_weights, diversification_ratio,
    min_variance_weights, risk_parity_weights,
)


COV = [[0.04, 0.01, 0.0], [0.01, 0.09, 0.02], [0.0, 0.02, 0.16]]


def test_weights_sum_to_one():
    w = max_diversification_weights(COV)
    assert sum(w) == pytest.approx(1.0, abs=1e-9)


def test_maximizes_diversification_ratio():
    w = max_diversification_weights(COV)
    dmax = diversification_ratio(w, COV)
    for other in (min_variance_weights(COV), risk_parity_weights(COV), [1/3]*3):
        assert diversification_ratio(other, COV) <= dmax + 1e-9


def test_ratio_at_least_one():
    w = max_diversification_weights(COV)
    assert diversification_ratio(w, COV) >= 1.0


def test_ratio_near_one_for_high_correlation():
    # Near-perfectly-correlated assets can't diversify: ratio ~ 1.
    hc = [[0.04, 0.0399], [0.0399, 0.04]]
    assert diversification_ratio([0.5, 0.5], hc) == pytest.approx(1.0, abs=1e-2)


def test_single_asset_ratio_is_one():
    assert diversification_ratio([1.0], [[0.04]]) == pytest.approx(1.0, abs=1e-12)


def test_validation():
    with pytest.raises(ValueError):
        max_diversification_weights([[1.0, 1.0], [1.0, 1.0]])   # singular
