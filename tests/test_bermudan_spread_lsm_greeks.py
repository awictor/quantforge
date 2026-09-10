"""Greeks of an American spread option by CRN bumps (bermudan_spread_lsm_greeks)."""

import pytest

from quantforge import (
    bermudan_spread_lsm_greeks,
    spread_greeks,
    OptionType,
)


S1, S2, K, T, R = 100.0, 95.0, 5.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.25, 0.3, 0.4


@pytest.mark.slow
def test_deltas_near_european_kirk():
    # Without dividends the American spread call ~ European, so its spot deltas
    # are close to the Kirk spread_greeks deltas.
    g = bermudan_spread_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                   n_steps=40, n_paths=120_000, seed=1)
    eg = spread_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO)
    assert g["delta1"] == pytest.approx(eg["delta1"], abs=0.03)
    assert g["delta2"] == pytest.approx(eg["delta2"], abs=0.03)


def test_long_leg_positive_short_leg_negative():
    g = bermudan_spread_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                   n_steps=12, n_paths=5_000, seed=2)
    assert g["delta1"] > 0.0   # long the spread in S1
    assert g["delta2"] < 0.0   # short the spread in S2


def test_put_delta_signs_flip():
    g = bermudan_spread_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                   option_type=OptionType.PUT, n_steps=12,
                                   n_paths=5_000, seed=3)
    assert g["delta1"] < 0.0
    assert g["delta2"] > 0.0


def test_reproducible():
    kw = dict(n_steps=6, n_paths=800, seed=99)
    a = bermudan_spread_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, **kw)
    b = bermudan_spread_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, **kw)
    assert a["delta1"] == b["delta1"]
    assert a["price"] == b["price"]
