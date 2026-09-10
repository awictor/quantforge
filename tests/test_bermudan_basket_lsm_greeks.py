"""Greeks of an American basket option by CRN bumps (bermudan_basket_lsm_greeks)."""

import pytest

from quantforge import (
    bermudan_basket_lsm_greeks,
    basket_greeks,
    OptionType,
)


S1, S2, K, T, R = 100.0, 90.0, 95.0, 1.0, 0.03
W1, W2 = 0.6, 0.4
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


@pytest.mark.slow
def test_deltas_near_european():
    # Without dividends the American basket call ~ European, so its spot deltas
    # are close to the moment-matched basket_greeks deltas.
    g = bermudan_basket_lsm_greeks(S1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO,
                                   n_steps=40, n_paths=120_000, seed=1)
    eg = basket_greeks((S1, S2), (W1, W2), K, T, R, (SIG1, SIG2), RHO)
    assert g["delta1"] == pytest.approx(eg["delta1"], abs=0.03)
    assert g["delta2"] == pytest.approx(eg["delta2"], abs=0.03)


def test_call_deltas_positive():
    g = bermudan_basket_lsm_greeks(S1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO,
                                   n_steps=12, n_paths=3_000, seed=2)
    assert g["delta1"] > 0.0
    assert g["delta2"] > 0.0


def test_higher_weight_larger_delta():
    # Asset 1 carries the larger basket weight, so the basket is more sensitive
    # to it per unit spot -> larger delta1 than delta2 (spots are comparable).
    g = bermudan_basket_lsm_greeks(S1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO,
                                   n_steps=15, n_paths=6_000, seed=3)
    assert g["delta1"] > g["delta2"]


def test_put_deltas_negative():
    g = bermudan_basket_lsm_greeks(S1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO,
                                   option_type=OptionType.PUT, n_steps=12,
                                   n_paths=3_000, seed=4)
    assert g["delta1"] < 0.0
    assert g["delta2"] < 0.0


def test_reproducible():
    kw = dict(n_steps=15, n_paths=4_000, seed=99)
    a = bermudan_basket_lsm_greeks(S1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO, **kw)
    b = bermudan_basket_lsm_greeks(S1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO, **kw)
    assert a["delta1"] == b["delta1"]
    assert a["price"] == b["price"]
