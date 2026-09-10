"""American basket option by LSM (bermudan_basket_lsm)."""

import pytest

from quantforge import bermudan_basket_lsm, basket_option, OptionType


S1, S2, K, T, R = 100.0, 90.0, 95.0, 1.0, 0.03
W1, W2 = 0.6, 0.4
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


@pytest.mark.slow
def test_no_dividend_call_close_to_european():
    eu = basket_option((S1, S2), (W1, W2), K, T, R, (SIG1, SIG2), RHO,
                       option_type=OptionType.CALL)
    am = bermudan_basket_lsm(S1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO,
                             n_steps=50, n_paths=60_000, seed=1)
    assert am == pytest.approx(eu, abs=0.1)


@pytest.mark.slow
def test_dividend_call_has_premium():
    q = 0.06
    eu = basket_option((S1, S2), (W1, W2), K, T, R, (SIG1, SIG2), RHO,
                       q=(q, q), option_type=OptionType.CALL)
    am = bermudan_basket_lsm(S1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO, q, q,
                             n_steps=50, n_paths=60_000, seed=2)
    assert am > eu + 0.1


@pytest.mark.slow
def test_put_has_early_exercise_premium():
    eu = basket_option((S1, S2), (W1, W2), K, T, R, (SIG1, SIG2), RHO,
                       option_type=OptionType.PUT)
    am = bermudan_basket_lsm(S1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO,
                             option_type=OptionType.PUT, n_steps=50,
                             n_paths=60_000, seed=3)
    assert am > eu + 0.05


def test_price_positive():
    am = bermudan_basket_lsm(S1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO,
                             n_steps=15, n_paths=6_000, seed=4)
    assert am > 0.0


def test_reproducible():
    kw = dict(n_steps=15, n_paths=4_000, seed=99)
    a = bermudan_basket_lsm(S1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO, **kw)
    b = bermudan_basket_lsm(S1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO, **kw)
    assert a == b


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        bermudan_basket_lsm(S1, S2, W1, W2, K, T, R, SIG1, SIG2, 1.5)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        bermudan_basket_lsm(-1, S2, W1, W2, K, T, R, SIG1, SIG2, RHO)
