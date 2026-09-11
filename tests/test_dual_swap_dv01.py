"""Risk of a dual-curve swap (dual_swap_dv01)."""

import math

import pytest

from quantforge import dual_swap_dv01, dual_swap_value, DiscountCurve


TS = [1, 2, 3, 4, 5]
OIS = DiscountCurve(TS, [math.exp(-0.025 * t) for t in TS])
PROJ = DiscountCurve(TS, [math.exp(-0.03 * t) for t in TS])
PAY = [1, 2, 3, 4, 5]
K = 0.03


def test_pv01_equals_minus_dv_dfixed_for_payer():
    g = dual_swap_dv01(OIS, PROJ, PAY, K, payer=True)
    h = 1e-6
    fd = (dual_swap_value(OIS, PROJ, PAY, K + h, True)
          - dual_swap_value(OIS, PROJ, PAY, K - h, True)) / (2 * h)
    assert g["pv01"] == pytest.approx(-fd, abs=1e-6)


def test_dv01_additive_across_curves():
    # To first order the total parallel DV01 is the sum of the OIS-only and
    # projection-only shifts.
    g = dual_swap_dv01(OIS, PROJ, PAY, K, payer=True)
    assert g["dv01"] == pytest.approx(g["ois_dv01"] + g["proj_dv01"], abs=1e-6)


def test_payer_dv01_negative_on_rate_drop():
    # dv01 is the value change for a 1bp *drop*; a payer loses when rates fall.
    g = dual_swap_dv01(OIS, PROJ, PAY, K, payer=True)
    assert g["dv01"] < 0.0


def test_receiver_dv01_positive_on_rate_drop():
    g = dual_swap_dv01(OIS, PROJ, PAY, K, payer=False)
    assert g["dv01"] > 0.0


def test_pv01_positive():
    g = dual_swap_dv01(OIS, PROJ, PAY, K, payer=True)
    assert g["pv01"] > 0.0


def test_value_field_matches_swap_value():
    g = dual_swap_dv01(OIS, PROJ, PAY, K, payer=True)
    assert g["value"] == pytest.approx(
        dual_swap_value(OIS, PROJ, PAY, K, True), abs=1e-12)
