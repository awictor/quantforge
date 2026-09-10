"""Greeks of a compound option (Geske) by FD (compound_option_greeks)."""

import pytest

from quantforge import compound_option_greeks, compound_option


S, K1, K2, T1, T2, R, SIG = 100.0, 5.0, 100.0, 0.5, 1.0, 0.05, 0.2


def test_delta_matches_finite_difference():
    g = compound_option_greeks(S, K1, K2, T1, T2, R, SIG, "call-on-call")
    h = 0.01
    fd = (compound_option(S + h, K1, K2, T1, T2, R, SIG, "call-on-call")
          - compound_option(S - h, K1, K2, T1, T2, R, SIG, "call-on-call")) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-5)


def test_call_on_call_greek_signs():
    g = compound_option_greeks(S, K1, K2, T1, T2, R, SIG, "call-on-call")
    assert g["delta"] > 0.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_put_on_call_delta_negative():
    # A put on a call loses value as the underlying call gains -> negative delta.
    g = compound_option_greeks(S, K1, K2, T1, T2, R, SIG, "put-on-call")
    assert g["delta"] < 0.0


def test_price_field_matches_compound_option():
    g = compound_option_greeks(S, K1, K2, T1, T2, R, SIG, "call-on-call")
    assert g["price"] == pytest.approx(
        compound_option(S, K1, K2, T1, T2, R, SIG, "call-on-call"), abs=1e-12)


def test_call_compounds_have_positive_vega():
    # A call *on* an option is long the compound optionality -> positive vega.
    # A put on an option can have negative vega (it is short that optionality),
    # so the sign is not universal across kinds.
    for kind in ("call-on-call", "call-on-put"):
        g = compound_option_greeks(S, K1, K2, T1, T2, R, SIG, kind)
        assert g["vega"] > 0.0, kind


def test_vega_defined_for_all_kinds():
    # All four kinds return a finite vega (positive for calls, can be negative
    # for puts on an option).
    for kind in ("call-on-call", "call-on-put", "put-on-call", "put-on-put"):
        g = compound_option_greeks(S, K1, K2, T1, T2, R, SIG, kind)
        assert g["vega"] == g["vega"]  # not NaN


def test_bad_kind_raises():
    with pytest.raises(ValueError):
        compound_option_greeks(S, K1, K2, T1, T2, R, SIG, "call-on-swap")
