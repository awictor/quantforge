"""Double-Heston (two-factor stochastic-vol) Fourier pricer."""

import math

import pytest

from quantforge import (
    OptionType,
    double_heston_price,
    double_heston_smile,
    heston_price,
)


HESTON_SETS = [
    (100, 100, 1.0, 0.03, 0.0, 0.04, 1.5, 0.04, 0.5, -0.7),
    (100, 120, 2.0, 0.05, 0.02, 0.09, 1.0, 0.06, 1.0, -0.5),
    (100, 90, 0.5, 0.04, 0.0, 0.04, 2.0, 0.04, 0.4, -0.6),
]


@pytest.mark.parametrize("S,K,t,r,q,v0,kappa,theta,xi,rho", HESTON_SETS)
def test_second_factor_off_recovers_heston(S, K, t, r, q, v0, kappa, theta, xi, rho):
    dh = double_heston_price(S, K, t, r, v0, kappa, theta, xi, rho,
                             0.0, 1.0, 0.0, 0.0, 0.0, OptionType.CALL, q=q)
    h = heston_price(S, K, t, r, v0, kappa, theta, xi, rho, OptionType.CALL, q=q)
    assert dh == pytest.approx(h, abs=1e-8)


def test_put_call_parity():
    S, K, t, r, q = 100, 105, 1.0, 0.04, 0.01
    f = (1.5, 0.04, 0.5, -0.7, 0.3, 0.02, 0.3, -0.4)
    c = double_heston_price(S, K, t, r, 0.03, 1.5, 0.04, 0.5, -0.7,
                            0.02, 0.3, 0.02, 0.3, -0.4, OptionType.CALL, q=q)
    p = double_heston_price(S, K, t, r, 0.03, 1.5, 0.04, 0.5, -0.7,
                            0.02, 0.3, 0.02, 0.3, -0.4, OptionType.PUT, q=q)
    rhs = S * math.exp(-q * t) - K * math.exp(-r * t)
    assert (c - p) == pytest.approx(rhs, abs=1e-8)


def test_two_factors_worth_more_than_one():
    # Adding a second variance factor adds variance, so an ATM call is dearer.
    S, K, t, r = 100, 100, 1.0, 0.03
    one = double_heston_price(S, K, t, r, 0.04, 2.0, 0.04, 0.4, -0.7,
                              0.0, 1.0, 0.0, 0.0, 0.0, OptionType.CALL)
    two = double_heston_price(S, K, t, r, 0.04, 2.0, 0.04, 0.4, -0.7,
                              0.03, 0.5, 0.04, 0.3, -0.5, OptionType.CALL)
    assert two > one


def test_downward_skew_from_negative_correlations():
    sm = double_heston_smile(100, [85, 92, 100, 108, 116], 1.0, 0.03,
                             0.04, 2.0, 0.04, 0.4, -0.7,
                             0.04, 0.5, 0.05, 0.3, -0.5)
    vols = [iv for _, iv in sm]
    assert vols[0] > vols[-1]


def test_intrinsic_at_expiry():
    v = double_heston_price(100, 90, 0.0, 0.03, 0.04, 1.5, 0.04, 0.5, -0.7,
                            0.02, 0.3, 0.02, 0.3, -0.4, OptionType.CALL)
    assert v == pytest.approx(10.0)


@pytest.mark.slow
def test_matches_qe_monte_carlo():
    import random
    from quantforge.heston_mc import _norm_ppf

    f1 = (0.04, 2.0, 0.04, 0.4, -0.7)
    f2 = (0.04, 0.5, 0.05, 0.3, -0.5)
    S, K, t, r, q = 100, 100, 1.0, 0.03, 0.0
    cf = double_heston_price(S, K, t, r, *f1, *f2, OptionType.CALL, q=q)

    ns, npths, seed = 100, 120000, 7
    dt = t / ns
    disc = math.exp(-r * t)
    rng = random.Random(seed)
    sm = []
    for _ in range(npths // 2):
        z = [[rng.gauss(0, 1), rng.gauss(0, 1)] for _ in range(ns)]
        u = [[rng.random(), rng.random()] for _ in range(ns)]
        for sgn in (1, -1):
            x = math.log(S)
            vs = [f1[0], f2[0]]
            for i in range(ns):
                incr = (r - q) * dt
                for fi, f in enumerate((f1, f2)):
                    v0, kappa, theta, xi, rho = f
                    ekt = math.exp(-kappa * dt)
                    m = theta + (vs[fi] - theta) * ekt
                    s2 = (vs[fi] * xi * xi * ekt / kappa) * (1 - ekt) \
                        + (theta * xi * xi / (2 * kappa)) * (1 - ekt) ** 2
                    uz = u[i][fi] if sgn == 1 else 1 - u[i][fi]
                    if m <= 0:
                        vn = 0.0
                    else:
                        psi = s2 / (m * m)
                        if psi <= 1.5:
                            inv = 2 / psi
                            b2 = inv - 1 + math.sqrt(inv) * math.sqrt(max(inv - 1, 0))
                            aa = m / (1 + b2)
                            vn = aa * (math.sqrt(b2) + _norm_ppf(uz)) ** 2
                        else:
                            pp = (psi - 1) / (psi + 1)
                            beta = (1 - pp) / m
                            vn = 0.0 if uz <= pp else math.log((1 - pp) / (1 - uz)) / beta
                    K0 = -rho * kappa * theta * dt / xi
                    K1 = 0.5 * dt * (kappa * rho / xi - 0.5) - rho / xi
                    K2 = 0.5 * dt * (kappa * rho / xi - 0.5) + rho / xi
                    K3 = 0.5 * dt * (1 - rho * rho)
                    K4 = 0.5 * dt * (1 - rho * rho)
                    vol = math.sqrt(max(K3 * vs[fi] + K4 * vn, 0))
                    incr += K0 + K1 * vs[fi] + K2 * vn + vol * (sgn * z[i][fi])
                    vs[fi] = vn
                x += incr
            sm.append(disc * max(math.exp(x) - K, 0))
    mean = sum(sm) / len(sm)
    var = sum((a - mean) ** 2 for a in sm) / (len(sm) - 1)
    se = math.sqrt(var / len(sm))
    assert abs(cf - mean) < 4.0 * se + 0.02
