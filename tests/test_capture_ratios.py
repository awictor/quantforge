"""Up/down capture ratios and downside beta (perfmetrics module)."""

import random

import pytest

from quantforge import up_capture, down_capture, downside_beta


def _market(n=400, seed=3):
    rng = random.Random(seed)
    return [rng.gauss(0.0, 0.01) for _ in range(n)]


def test_self_capture_is_one():
    m = _market()
    assert up_capture(m, m) == pytest.approx(1.0, abs=1e-9)
    assert down_capture(m, m) == pytest.approx(1.0, abs=1e-9)
    assert downside_beta(m, m) == pytest.approx(1.0, abs=1e-9)


def test_leveraged_captures_above_one():
    m = _market()
    lev = [1.5 * x for x in m]
    assert up_capture(lev, m) == pytest.approx(1.5, abs=0.02)
    assert down_capture(lev, m) == pytest.approx(1.5, abs=0.02)
    assert downside_beta(lev, m) == pytest.approx(1.5, abs=1e-6)


def test_defensive_captures_below_one():
    m = _market()
    defn = [0.5 * x for x in m]
    assert up_capture(defn, m) < 1.0
    assert down_capture(defn, m) < 1.0


def test_downside_beta_uses_only_down_periods():
    m = _market()
    # An asset that mirrors the market only on down days.
    a = [x if x < 0 else 0.0 for x in m]
    assert downside_beta(a, m) == pytest.approx(1.0, abs=1e-9)


def test_validation():
    m = _market(n=50)
    with pytest.raises(ValueError):
        up_capture(m, m[:10])                 # length mismatch
    with pytest.raises(ValueError):
        up_capture([0.01, 0.02], [-0.01, -0.02])  # no up periods
    with pytest.raises(ValueError):
        downside_beta([0.01, 0.02], [0.01, 0.02])  # no down periods
