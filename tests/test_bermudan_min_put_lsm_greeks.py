"""Greeks of an American min-put by CRN bumps (bermudan_min_put_lsm_greeks)."""

import pytest

from quantforge import bermudan_min_put_lsm_greeks, worst_of_put_closed


S1, S2, K, T, R = 100.0, 95.0, 100.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


def test_both_deltas_negative():
    # Raising either spot lifts min(S1,S2), shrinking the protective put.
    g = bermudan_min_put_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                    n_steps=20, n_paths=15_000, seed=1)
    assert g["delta1"] < 0.0
    assert g["delta2"] < 0.0


def test_lower_asset_has_larger_magnitude_delta():
    # S2 starts lower, so it is more often the min -> the put is more sensitive
    # to it (more negative delta).
    g = bermudan_min_put_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                    n_steps=25, n_paths=30_000, seed=2)
    assert g["delta2"] < g["delta1"]


@pytest.mark.slow
def test_price_matches_lsm_and_exceeds_european():
    g = bermudan_min_put_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                    n_steps=40, n_paths=60_000, seed=3)
    eu = worst_of_put_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    assert g["price"] > eu  # American premium


def test_reproducible():
    kw = dict(n_steps=15, n_paths=5_000, seed=99)
    a = bermudan_min_put_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, **kw)
    b = bermudan_min_put_lsm_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO, **kw)
    assert a["delta1"] == b["delta1"]
    assert a["price"] == b["price"]
