"""Tests for the American early-exercise premium decomposition."""

import pytest

from quantforge import (
    early_exercise_premium, bjerksund_stensland, call_price, put_price, OptionType,
)


def test_no_dividend_call_premium_is_zero():
    # b = r: American call == European call, so the premium is zero.
    d = early_exercise_premium(100, 100, 1.0, 0.05, 0.2, OptionType.CALL)
    assert d["premium"] == pytest.approx(0.0, abs=1e-6)
    assert d["american"] == pytest.approx(d["european"], abs=1e-6)


def test_decomposition_adds_up():
    d = early_exercise_premium(90, 100, 1.0, 0.05, 0.3, OptionType.PUT)
    assert d["american"] == pytest.approx(d["european"] + d["premium"], abs=1e-9)


def test_itm_put_has_positive_premium():
    d = early_exercise_premium(90, 100, 1.0, 0.05, 0.3, OptionType.PUT)
    assert d["premium"] > 0
    assert d["european"] == pytest.approx(put_price(90, 100, 1.0, 0.05, 0.3), abs=1e-9)


def test_dividend_call_has_positive_premium():
    d = early_exercise_premium(100, 100, 1.0, 0.05, 0.25, OptionType.CALL, b=-0.05)
    assert d["premium"] > 0


def test_american_matches_bjerksund_stensland():
    d = early_exercise_premium(95, 100, 0.5, 0.06, 0.3, OptionType.PUT, b=0.03)
    assert d["american"] == pytest.approx(
        bjerksund_stensland(95, 100, 0.5, 0.06, 0.3, OptionType.PUT, b=0.03), abs=1e-9)


def test_premium_grows_with_dividend():
    small = early_exercise_premium(100, 100, 1.0, 0.05, 0.25, OptionType.CALL, b=0.0)
    large = early_exercise_premium(100, 100, 1.0, 0.05, 0.25, OptionType.CALL, b=-0.10)
    assert large["premium"] > small["premium"]
