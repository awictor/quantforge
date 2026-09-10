"""SABR implied-density arbitrage detection and repair."""

import pytest

from quantforge import (
    SABRParams,
    sabr_density,
    sabr_butterfly_arbitrage,
    sabr_is_arbitrage_free,
    sabr_repair_butterfly,
)


F, T = 0.04, 1.0
BENIGN = (0.02, 0.5, -0.3, 0.4)     # alpha, beta, rho, nu
ARB = (0.02, 0.5, -0.9, 2.5)        # steep wings -> negative density


def test_density_integrates_to_one():
    ks = [F * (0.1 + 0.005 * i) for i in range(600)]
    dK = ks[1] - ks[0]
    integ = sum(sabr_density(F, K, T, *BENIGN) for K in ks) * dK
    assert integ == pytest.approx(1.0, abs=1e-2)


def test_benign_smile_is_arbitrage_free():
    assert sabr_is_arbitrage_free(F, T, *BENIGN)
    assert sabr_butterfly_arbitrage(F, T, *BENIGN) == []


def test_atm_density_positive():
    assert sabr_density(F, F, T, *BENIGN) > 0.0


def test_extreme_vol_of_vol_flags_arbitrage():
    bad = sabr_butterfly_arbitrage(F, T, *ARB)
    assert len(bad) > 0
    assert not sabr_is_arbitrage_free(F, T, *ARB)


def test_repair_restores_arbitrage_freedom():
    params = SABRParams(alpha=ARB[0], beta=ARB[1], rho=ARB[2], nu=ARB[3])
    fixed = sabr_repair_butterfly(F, T, params)
    assert fixed.nu < params.nu
    assert sabr_is_arbitrage_free(F, T, fixed.alpha, fixed.beta, fixed.rho,
                                  fixed.nu)
    # Level, skew and elasticity are preserved.
    assert (fixed.alpha, fixed.beta, fixed.rho) == (params.alpha, params.beta,
                                                    params.rho)


def test_repair_leaves_clean_smile_unchanged():
    params = SABRParams(alpha=BENIGN[0], beta=BENIGN[1], rho=BENIGN[2],
                        nu=BENIGN[3])
    fixed = sabr_repair_butterfly(F, T, params)
    assert fixed.nu == params.nu
