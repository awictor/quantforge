"""American call on the max of two assets by LSM (bermudan_max_call_lsm)."""

import pytest

from quantforge import bermudan_max_call_lsm, best_of_call_closed


S1, S2, K, T, R = 100.0, 95.0, 100.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


@pytest.mark.slow
def test_no_dividend_close_to_european():
    # Without dividends the max-call has little early-exercise value, so the LSM
    # American price is close to the European Stulz closed form.
    eu = best_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    am = bermudan_max_call_lsm(S1, S2, K, T, R, SIG1, SIG2, RHO,
                               n_steps=50, n_paths=60_000, seed=1)
    assert am == pytest.approx(eu, rel=0.02)


@pytest.mark.slow
def test_dividends_give_early_exercise_premium():
    q = 0.06
    eu = best_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, RHO, q, q)
    am = bermudan_max_call_lsm(S1, S2, K, T, R, SIG1, SIG2, RHO, q, q,
                               n_steps=50, n_paths=60_000, seed=2)
    # American must exceed European by a positive early-exercise premium.
    assert am > eu + 0.05


def test_price_positive():
    am = bermudan_max_call_lsm(S1, S2, K, T, R, SIG1, SIG2, RHO,
                               n_steps=20, n_paths=8_000, seed=3)
    assert am > 0.0


def test_at_least_single_asset_american_lower_bound():
    # The max-call is worth at least a European vanilla on either asset.
    from quantforge import call_price
    am = bermudan_max_call_lsm(S1, S2, K, T, R, SIG1, SIG2, RHO,
                               n_steps=30, n_paths=20_000, seed=4)
    assert am > call_price(S1, K, T, R, SIG1) - 0.5


def test_reproducible():
    kw = dict(n_steps=20, n_paths=5_000, seed=99)
    a = bermudan_max_call_lsm(S1, S2, K, T, R, SIG1, SIG2, RHO, **kw)
    b = bermudan_max_call_lsm(S1, S2, K, T, R, SIG1, SIG2, RHO, **kw)
    assert a == b


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        bermudan_max_call_lsm(S1, S2, K, T, R, SIG1, SIG2, 1.5)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        bermudan_max_call_lsm(-1, S2, K, T, R, SIG1, SIG2, RHO)
