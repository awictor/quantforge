"""Kou double-exponential jump-diffusion characteristic-function pricer."""

import math

import pytest

from quantforge import OptionType, kou_price, kou_smile, call_price


def test_zero_intensity_recovers_black_scholes():
    for S, K, t, r, sig, q in [(100, 100, 1.0, 0.05, 0.2, 0.0),
                               (100, 110, 0.5, 0.03, 0.3, 0.02)]:
        k0 = kou_price(S, K, t, r, sig, 0.0, 0.5, 10.0, 10.0, OptionType.CALL, q=q)
        bs = call_price(S, K, t, r, sig, b=r - q)
        assert k0 == pytest.approx(bs, abs=1e-10)


def test_put_call_parity():
    S, K, t, r, q = 100, 105, 1.0, 0.04, 0.01
    args = (0.2, 3.0, 0.4, 10.0, 5.0)
    c = kou_price(S, K, t, r, *args, OptionType.CALL, q=q)
    p = kou_price(S, K, t, r, *args, OptionType.PUT, q=q)
    rhs = S * math.exp(-q * t) - K * math.exp(-r * t)
    assert (c - p) == pytest.approx(rhs, abs=1e-8)


def test_jumps_raise_price_vs_black_scholes():
    S, K, t, r, sig = 100, 100, 0.5, 0.03, 0.15
    bs = call_price(S, K, t, r, sig, b=r)
    kou = kou_price(S, K, t, r, sig, 3.0, 0.4, 10.0, 5.0, OptionType.CALL)
    assert kou > bs


def test_downward_skew_from_asymmetric_jumps():
    # More/heavier down jumps (eta2 < eta1, p < 1/2) => downside vol richer.
    S, t, r = 100.0, 0.25, 0.03
    strikes = [85, 92, 100, 108, 116]
    sm = kou_smile(S, strikes, t, r, 0.15, 4.0, 0.3, 15.0, 5.0)
    vols = [iv for _, iv in sm]
    assert vols[0] > vols[-1]


def test_intrinsic_at_expiry():
    assert kou_price(100, 90, 0.0, 0.03, 0.2, 3.0, 0.4, 10.0, 5.0,
                     OptionType.CALL) == pytest.approx(10.0)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        kou_price(100, 100, 1.0, 0.05, 0.2, 1.0, 0.5, 0.8, 5.0)  # eta1 <= 1
    with pytest.raises(ValueError):
        kou_price(100, 100, 1.0, 0.05, 0.2, 1.0, 1.5, 10.0, 5.0)  # p out of [0,1]
    with pytest.raises(ValueError):
        kou_price(100, 100, 1.0, 0.05, 0.2, -1.0, 0.5, 10.0, 5.0)  # lam < 0


@pytest.mark.slow
@pytest.mark.parametrize("S,K,t,r,sigma,lam,p,eta1,eta2,q", [
    (100, 100, 0.5, 0.05, 0.15, 3.0, 0.4, 10.0, 5.0, 0.0),
    (100, 90, 1.0, 0.03, 0.2, 2.0, 0.3, 8.0, 4.0, 0.02),
    (100, 110, 0.25, 0.04, 0.1, 5.0, 0.5, 15.0, 10.0, 0.0),
])
def test_char_function_matches_monte_carlo(S, K, t, r, sigma, lam, p,
                                           eta1, eta2, q):
    import random

    zeta = p * eta1 / (eta1 - 1.0) + (1.0 - p) * eta2 / (eta2 + 1.0) - 1.0
    drift = (r - q - 0.5 * sigma * sigma - lam * zeta) * t
    vol = sigma * math.sqrt(t)
    disc = math.exp(-r * t)
    rng = random.Random(7)
    sm = []
    for _ in range(500_000):
        x = drift + vol * rng.gauss(0, 1)
        cnt, term, cum, uu = 0, math.exp(-lam * t), math.exp(-lam * t), rng.random()
        while uu > cum and cnt < 200:
            cnt += 1
            term *= lam * t / cnt
            cum += term
        for _j in range(cnt):
            x += rng.expovariate(eta1) if rng.random() < p else -rng.expovariate(eta2)
        sm.append(disc * max(S * math.exp(x) - K, 0.0))
    mean = sum(sm) / len(sm)
    var = sum((a - mean) ** 2 for a in sm) / (len(sm) - 1)
    se = math.sqrt(var / len(sm))

    cf = kou_price(S, K, t, r, sigma, lam, p, eta1, eta2, OptionType.CALL, q=q)
    assert abs(cf - mean) < 4.0 * se + 0.01
