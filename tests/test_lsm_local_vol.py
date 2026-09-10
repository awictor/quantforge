"""Longstaff-Schwartz American pricing on a local-volatility surface."""

import pytest

from quantforge import (
    OptionType,
    bermudan_lsm_local_vol,
    bermudan_lsm,
    american_price as crr,
    crank_nicolson_price,
)


S, K, T, R, SIGMA = 100.0, 100.0, 1.0, 0.05, 0.2


@pytest.mark.slow
def test_flat_local_vol_matches_constant_lsm_and_tree():
    lv = bermudan_lsm_local_vol(S, K, T, R, lambda s, tau: SIGMA,
                                OptionType.PUT, n_steps=50, n_paths=40_000,
                                seed=1)
    const = bermudan_lsm(S, K, T, R, SIGMA, OptionType.PUT, n_steps=50,
                         n_paths=40_000, seed=1)
    tree = crr(S, K, T, R, SIGMA, OptionType.PUT, b=R, steps=2000)
    assert lv == pytest.approx(const, abs=1e-6)
    assert lv == pytest.approx(tree, abs=5e-2)


@pytest.mark.slow
def test_skewed_local_vol_matches_pde_american():
    def lv(s, tau):
        return max(0.1, 0.25 - 0.001 * (s - 100))
    am_lsm = bermudan_lsm_local_vol(S, K, T, R, lv, OptionType.PUT,
                                    n_steps=80, n_paths=60_000, seed=2)
    am_pde = crank_nicolson_price(S, K, T, R, option_type=OptionType.PUT,
                                  american=True,
                                  local_vol_fn=lambda s, tt: max(0.1, 0.25 - 0.001 * (s - 100)),
                                  n_space=300, n_time=300)
    assert am_lsm == pytest.approx(am_pde, abs=0.1)


def test_american_put_above_european_ish():
    # American LSM >= the terminal-only (European) LSM proxy at n_steps=1.
    am = bermudan_lsm_local_vol(S, K, T, R, lambda s, tau: SIGMA,
                                OptionType.PUT, n_steps=25, n_paths=8_000,
                                seed=3)
    eu = bermudan_lsm_local_vol(S, K, T, R, lambda s, tau: SIGMA,
                                OptionType.PUT, n_steps=1, n_paths=8_000,
                                seed=3)
    assert am >= eu - 0.1


def test_bad_params_raise():
    with pytest.raises(ValueError):
        bermudan_lsm_local_vol(-1, K, T, R, lambda s, tau: SIGMA)
    with pytest.raises(ValueError):
        bermudan_lsm_local_vol(S, K, T, R, lambda s, tau: SIGMA, n_steps=0)
