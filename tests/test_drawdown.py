"""Drawdown duration and recovery analytics."""

from quantforge import drawdown_analytics, max_drawdown


def test_depth_matches_max_drawdown():
    rets = [0.1, 0.1, -0.2, -0.125, 0.180, 0.30]
    a = drawdown_analytics(rets)
    assert abs(a["max_drawdown_depth"] - max_drawdown(rets)) < 1e-12


def test_peak_before_trough_before_recovery():
    rets = [0.1, 0.1, -0.2, -0.125, 0.180, 0.30]
    a = drawdown_analytics(rets)
    assert a["peak_index"] < a["trough_index"] <= a["recovery_index"]
    assert a["time_to_recovery"] == a["recovery_index"] - a["trough_index"]


def test_monotone_up_has_no_drawdown():
    a = drawdown_analytics([0.05] * 10)
    assert a["max_drawdown_depth"] == 0.0
    assert a["longest_underwater"] == 0
    assert a["time_to_recovery"] == 0


def test_never_recovers():
    a = drawdown_analytics([0.1, 0.1, -0.5, -0.1])
    assert a["recovery_index"] is None
    assert a["time_to_recovery"] is None


def test_longest_underwater_counts_periods():
    rets = [0.1, -0.05, 0.06, 0.1, -0.02, -0.02, -0.02, -0.02, 0.20]
    a = drawdown_analytics(rets)
    assert a["longest_underwater"] == 4


def test_empty_series():
    a = drawdown_analytics([])
    assert a["max_drawdown_depth"] == 0.0
    assert a["time_to_recovery"] == 0
