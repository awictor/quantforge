"""Variance-swap fair strike replicated from a volatility smile."""

import math

import pytest

from quantforge import variance_swap_from_smile, variance_swap_strike


S0, T, R = 100.0, 1.0, 0.03


@pytest.mark.parametrize("sig", [0.15, 0.2, 0.3])
def test_flat_smile_gives_variance(sig):
    # Model-free: a flat smile's variance-swap strike equals sigma^2.
    kv = variance_swap_from_smile(S0, T, R, lambda K: sig,
                                  n_strikes=401, width=8.0)
    assert kv == pytest.approx(sig * sig, abs=3e-3)


def test_converges_to_flat_variance_with_wider_strip():
    sig = 0.2
    e_narrow = abs(variance_swap_from_smile(S0, T, R, lambda K: sig,
                                            n_strikes=201, width=5.0) - 0.04)
    e_wide = abs(variance_swap_from_smile(S0, T, R, lambda K: sig,
                                          n_strikes=801, width=10.0) - 0.04)
    assert e_wide < e_narrow


def test_skew_adds_variance_over_atm():
    # A downward skew lifts the replicated variance above the ATM level (the
    # wings contribute convexity).
    atm = 0.2

    def smile(K):
        # Bounded downward skew (kept positive in the far wings).
        return max(0.05, atm + 0.1 * math.log(S0 / K))

    kv = variance_swap_from_smile(S0, T, R, smile, n_strikes=201, width=4.0)
    # Compare to a flat-ATM smile on the *same* grid (both share the strip's
    # truncation), so the excess is purely the skew's convexity contribution.
    flat = variance_swap_from_smile(S0, T, R, lambda K: atm,
                                    n_strikes=201, width=4.0)
    assert kv > flat


def test_matches_manual_chain():
    # variance_swap_from_smile just builds a flat-vol chain and delegates; a
    # coarse manual chain at the same vol should give a close strike.
    from quantforge.bsm import call_price, put_price
    sig = 0.2
    F = S0 * math.exp(R * T)
    ks = [F * math.exp(x) for x in [i * 0.05 for i in range(-40, 41)]]
    pk = [K for K in ks if K < F]
    ck = [K for K in ks if K >= F]
    pp = [put_price(S0, K, T, R, sig, b=R) for K in pk]
    cp = [call_price(S0, K, T, R, sig, b=R) for K in ck]
    manual = variance_swap_strike(S0, T, R, pk, pp, ck, cp, split=F)
    smile_based = variance_swap_from_smile(S0, T, R, lambda K: sig,
                                           n_strikes=81, width=4.0)
    assert smile_based == pytest.approx(manual, abs=5e-3)


def test_bad_tenor_raises():
    with pytest.raises(ValueError):
        variance_swap_from_smile(S0, 0.0, R, lambda K: 0.2)
