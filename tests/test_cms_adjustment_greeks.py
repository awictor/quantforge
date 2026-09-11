"""Sensitivities of the standard-model CMS convexity adjustment (cms_adjustment_greeks)."""

import pytest

from quantforge import cms_adjustment_greeks, cms_adjustment_standard


F, SIG, EXP, TEN = 0.03, 0.2, 5.0, 10.0


def test_d_sigma_matches_finite_difference():
    g = cms_adjustment_greeks(F, SIG, EXP, TEN)
    h = 1e-6
    fd = (cms_adjustment_standard(F, SIG + h, EXP, TEN)
          - cms_adjustment_standard(F, SIG - h, EXP, TEN)) / (2 * h)
    assert g["d_sigma"] == pytest.approx(fd, abs=1e-6)


def test_d_forward_matches_finite_difference():
    g = cms_adjustment_greeks(F, SIG, EXP, TEN)
    h = 1e-6
    fd = (cms_adjustment_standard(F + h, SIG, EXP, TEN)
          - cms_adjustment_standard(F - h, SIG, EXP, TEN)) / (2 * h)
    assert g["d_forward"] == pytest.approx(fd, abs=1e-5)


def test_convexity_grows_with_vol():
    # More vol -> more convexity -> larger adjustment, so d_sigma > 0.
    g = cms_adjustment_greeks(F, SIG, EXP, TEN)
    assert g["d_sigma"] > 0.0


def test_adjustment_field_matches_standard():
    g = cms_adjustment_greeks(F, SIG, EXP, TEN)
    assert g["adjustment"] == pytest.approx(
        cms_adjustment_standard(F, SIG, EXP, TEN), abs=1e-12)


def test_zero_vol_zero_sensitivity():
    g = cms_adjustment_greeks(F, 1e-9, EXP, TEN)
    assert g["adjustment"] == pytest.approx(0.0, abs=1e-12)
