"""Greeks of rainbow options by FD on the Stulz closed forms (rainbow_greeks)."""

import pytest

from quantforge import rainbow_greeks, OptionType
from quantforge.bsm import delta as bs_delta


S1, S2, K, T, R = 100.0, 95.0, 100.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


def test_call_delta_identity():
    # d/dSi of the Stulz identity C_max + C_min = c(S1) + c(S2): the best-of and
    # worst-of call deltas in each asset sum to the single-asset BS delta.
    gb = rainbow_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, "best", OptionType.CALL)
    gw = rainbow_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, "worst", OptionType.CALL)
    assert gb["delta1"] + gw["delta1"] == pytest.approx(
        bs_delta(S1, K, T, R, SIG1, OptionType.CALL), abs=1e-4)
    assert gb["delta2"] + gw["delta2"] == pytest.approx(
        bs_delta(S2, K, T, R, SIG2, OptionType.CALL), abs=1e-4)


def test_corr_vega_signs_and_offset():
    # Max-call gains as correlation falls (corr_vega < 0); min-call gains as it
    # rises (corr_vega > 0); their sum is zero (the vanilla legs have no rho).
    gb = rainbow_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, "best", OptionType.CALL)
    gw = rainbow_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, "worst", OptionType.CALL)
    assert gb["corr_vega"] < 0.0
    assert gw["corr_vega"] > 0.0
    assert gb["corr_vega"] + gw["corr_vega"] == pytest.approx(0.0, abs=1e-3)


def test_deltas_in_unit_interval():
    gb = rainbow_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, "best", OptionType.CALL)
    assert 0.0 < gb["delta1"] < 1.0
    assert 0.0 < gb["delta2"] < 1.0


def test_own_gamma_positive():
    gb = rainbow_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, "best", OptionType.CALL)
    assert gb["gamma1"] > 0.0
    assert gb["gamma2"] > 0.0


def test_put_delta_negative():
    gw = rainbow_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, "worst", OptionType.PUT)
    assert gw["delta1"] < 0.0
    assert gw["delta2"] < 0.0


def test_put_delta_identity():
    gb = rainbow_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, "best", OptionType.PUT)
    gw = rainbow_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, "worst", OptionType.PUT)
    assert gb["delta1"] + gw["delta1"] == pytest.approx(
        bs_delta(S1, K, T, R, SIG1, OptionType.PUT), abs=1e-4)


def test_bad_kind_raises():
    with pytest.raises(ValueError):
        rainbow_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, kind="middle")


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        rainbow_greeks(S1, S2, K, T, R, SIG1, SIG2, 1.5)
