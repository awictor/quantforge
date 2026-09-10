"""Pathwise vega of an arithmetic-average Asian vs bump and closed form."""

import pytest

from quantforge import (
    OptionType,
    asian_pathwise_vega,
    arithmetic_asian,
)
from quantforge.montecarlo import arithmetic_asian_mc


S, K, T, R = 100.0, 100.0, 1.0, 0.05


@pytest.mark.slow
def test_pathwise_vega_matches_mc_bump():
    sig, ns = 0.3, 50
    pw = asian_pathwise_vega(S, K, T, R, sig, OptionType.CALL,
                             n_steps=ns, n_paths=200_000, seed=1)
    h = 1e-3
    up = arithmetic_asian_mc(S, K, T, R, sig + h, OptionType.CALL,
                             n_steps=ns, n_paths=200_000, seed=7)
    dn = arithmetic_asian_mc(S, K, T, R, sig - h, OptionType.CALL,
                             n_steps=ns, n_paths=200_000, seed=7)
    bump = (up.price - dn.price) / (2 * h)
    assert pw.price == pytest.approx(bump, abs=3.0 * pw.std_error + 0.5)


@pytest.mark.slow
def test_pathwise_vega_near_turnbull_wakeman_fd():
    h = 1e-3
    for sig in (0.2, 0.3):
        pw = asian_pathwise_vega(S, K, T, R, sig, OptionType.CALL,
                                 n_steps=100, n_paths=150_000, seed=2)
        tw = (arithmetic_asian(S, K, T, R, sig + h, OptionType.CALL)
              - arithmetic_asian(S, K, T, R, sig - h, OptionType.CALL)) / (2 * h)
        assert pw.price == pytest.approx(tw, abs=0.6)


def test_call_and_put_vega_positive():
    c = asian_pathwise_vega(S, K, T, R, 0.2, OptionType.CALL,
                            n_steps=20, n_paths=15_000, seed=3)
    p = asian_pathwise_vega(S, K, T, R, 0.2, OptionType.PUT,
                            n_steps=20, n_paths=15_000, seed=4)
    assert c.price > 0.0 and p.price > 0.0


def test_bad_params_raise():
    with pytest.raises(ValueError):
        asian_pathwise_vega(S, K, T, R, 0.2, n_steps=0)
