"""American put on the min of two assets by LSM (bermudan_min_put_lsm)."""

import pytest

from quantforge import bermudan_min_put_lsm, worst_of_put_closed


S1, S2, K, T, R = 100.0, 95.0, 100.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


@pytest.mark.slow
def test_exceeds_european_by_early_exercise_premium():
    # Puts carry early-exercise value even without dividends, so the American
    # min-put sits above the European Stulz worst_of_put_closed.
    eu = worst_of_put_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    am = bermudan_min_put_lsm(S1, S2, K, T, R, SIG1, SIG2, RHO,
                              n_steps=50, n_paths=60_000, seed=1)
    assert am > eu + 0.1


def test_at_least_european():
    eu = worst_of_put_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    am = bermudan_min_put_lsm(S1, S2, K, T, R, SIG1, SIG2, RHO,
                              n_steps=20, n_paths=10_000, seed=2)
    assert am > eu - 0.2


def test_price_positive():
    am = bermudan_min_put_lsm(S1, S2, K, T, R, SIG1, SIG2, RHO,
                              n_steps=15, n_paths=6_000, seed=3)
    assert am > 0.0


def test_deep_itm_near_intrinsic():
    # Very low spots: the put is deep ITM, worth about the discounted-free
    # intrinsic K - min(S1,S2) (immediate exercise optimal).
    am = bermudan_min_put_lsm(40.0, 45.0, 100.0, T, R, SIG1, SIG2, RHO,
                              n_steps=20, n_paths=8_000, seed=4)
    assert am >= 100.0 - 40.0 - 1e-6  # >= intrinsic on the min (40)


def test_reproducible():
    kw = dict(n_steps=15, n_paths=4_000, seed=99)
    a = bermudan_min_put_lsm(S1, S2, K, T, R, SIG1, SIG2, RHO, **kw)
    b = bermudan_min_put_lsm(S1, S2, K, T, R, SIG1, SIG2, RHO, **kw)
    assert a == b


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        bermudan_min_put_lsm(S1, S2, K, T, R, SIG1, SIG2, 1.5)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        bermudan_min_put_lsm(-1, S2, K, T, R, SIG1, SIG2, RHO)
