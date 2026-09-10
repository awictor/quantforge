"""Greeks of a two-asset correlated digital (two_asset_digital_greeks)."""

import pytest

from quantforge import two_asset_digital_greeks


S1, S2, K1, K2, T, R = 100.0, 95.0, 105.0, 90.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


def test_both_above_spot_deltas_positive():
    # A both-above digital pays more the higher either spot -> positive deltas.
    g = two_asset_digital_greeks(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                 "above", "above")
    assert g["delta1"] > 0.0
    assert g["delta2"] > 0.0


def test_correlation_vega_sign_by_quadrant():
    # Both-above (or both-below) gains as correlation rises; a mixed
    # above/below digital loses.
    both = two_asset_digital_greeks(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                    "above", "above")
    mixed = two_asset_digital_greeks(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                     "above", "below")
    assert both["corr_vega"] > 0.0
    assert mixed["corr_vega"] < 0.0


def test_correlation_vega_sums_to_zero():
    # Total probability is independent of rho, so the four quadrant corr-vegas
    # cancel.
    total = sum(
        two_asset_digital_greeks(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO, c1, c2)
        ["corr_vega"]
        for c1 in ("above", "below") for c2 in ("above", "below"))
    assert total == pytest.approx(0.0, abs=1e-6)


def test_cash_scales_greeks():
    a = two_asset_digital_greeks(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                 cash=1.0)
    b = two_asset_digital_greeks(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                 cash=5.0)
    assert b["delta1"] == pytest.approx(5.0 * a["delta1"], rel=1e-9)
    assert b["corr_vega"] == pytest.approx(5.0 * a["corr_vega"], rel=1e-9)


def test_below_above_put_side_delta1_negative():
    # cond1 = "below": the digital pays when S1 < K1, so it falls as S1 rises.
    g = two_asset_digital_greeks(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                 "below", "above")
    assert g["delta1"] < 0.0


def test_bad_condition_raises():
    with pytest.raises(ValueError):
        two_asset_digital_greeks(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                 cond1="maybe")
