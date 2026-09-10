"""Greeks of an American max-call by CRN bumps (bermudan_max_call_lsm_greeks)."""

import pytest

from quantforge import (
    bermudan_max_call_lsm_greeks,
    rainbow_greeks,
    OptionType,
)


S1, S2, K, T, R = 100.0, 95.0, 100.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


@pytest.mark.slow
def test_deltas_near_european_no_dividend():
    # Without dividends the American max-call ~ European, so its spot deltas are
    # close to the exact rainbow_greeks deltas.
    g = bermudan_max_call_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                     n_steps=40, n_paths=120_000, seed=1)
    eg = rainbow_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, "best",
                        OptionType.CALL)
    assert g["delta1"] == pytest.approx(eg["delta1"], abs=0.03)
    assert g["delta2"] == pytest.approx(eg["delta2"], abs=0.03)


def test_deltas_in_unit_interval():
    g = bermudan_max_call_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                     n_steps=10, n_paths=2_500, seed=2)
    assert 0.0 < g["delta1"] < 1.0
    assert 0.0 < g["delta2"] < 1.0


def test_cross_gamma_negative():
    # The two spots are substitutes in a max payoff: a higher S2 lowers the
    # sensitivity to S1, so the cross-gamma is negative.
    g = bermudan_max_call_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                     n_steps=12, n_paths=5_000, seed=3)
    assert g["cross"] < 0.0


def test_reproducible():
    kw = dict(n_steps=10, n_paths=1_500, seed=99)
    a = bermudan_max_call_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, **kw)
    b = bermudan_max_call_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, **kw)
    assert a["delta1"] == b["delta1"]
    assert a["price"] == b["price"]
