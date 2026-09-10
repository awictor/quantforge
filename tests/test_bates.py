"""Bates (Heston + Merton jumps) characteristic-function pricer."""

import math

import pytest

from quantforge import (
    OptionType,
    bates_price,
    bates_smile,
    heston_price,
    call_price,
)


# (S, K, t, r, v0, kappa, theta, xi, rho, q)
HESTON_SETS = [
    (100, 100, 1.0, 0.03, 0.04, 1.5, 0.04, 0.5, -0.7, 0.0),
    (100, 120, 2.0, 0.05, 0.09, 1.0, 0.06, 1.0, -0.5, 0.02),
    (100, 90, 0.5, 0.04, 0.04, 2.0, 0.04, 0.4, -0.6, 0.0),
]


@pytest.mark.parametrize("S,K,t,r,v0,kappa,theta,xi,rho,q", HESTON_SETS)
def test_zero_intensity_recovers_heston(S, K, t, r, v0, kappa, theta, xi, rho, q):
    bates = bates_price(S, K, t, r, v0, kappa, theta, xi, rho,
                        0.0, 0.0, 0.0, OptionType.CALL, q=q)
    hes = heston_price(S, K, t, r, v0, kappa, theta, xi, rho,
                       OptionType.CALL, q=q)
    assert bates == pytest.approx(hes, abs=1e-8)


def test_put_call_parity():
    S, K, t, r, q = 100, 105, 1.0, 0.04, 0.01
    args = (0.04, 1.5, 0.04, 0.6, -0.6, 0.5, -0.1, 0.15)
    c = bates_price(S, K, t, r, *args, OptionType.CALL, q=q)
    p = bates_price(S, K, t, r, *args, OptionType.PUT, q=q)
    rhs = S * math.exp(-q * t) - K * math.exp(-r * t)
    assert (c - p) == pytest.approx(rhs, abs=1e-8)


def test_jumps_raise_price_vs_heston():
    # Adding jump risk (nonzero lambda) increases an ATM call's value.
    S, K, t, r = 100, 100, 0.5, 0.03
    hes = heston_price(S, K, t, r, 0.04, 1.5, 0.04, 0.5, -0.7, OptionType.CALL)
    bat = bates_price(S, K, t, r, 0.04, 1.5, 0.04, 0.5, -0.7,
                      2.0, -0.1, 0.2, OptionType.CALL)
    assert bat > hes


def test_negative_mean_jump_steepens_downward_skew():
    S, t, r = 100.0, 0.25, 0.03
    strikes = [85, 92, 100, 108, 116]
    sm = bates_smile(S, strikes, t, r, 0.04, 1.5, 0.04, 0.4, -0.5,
                     3.0, -0.15, 0.15)
    vols = [iv for _, iv in sm]
    assert vols[0] > vols[-1]  # downside vol richer


def test_intrinsic_at_expiry():
    assert bates_price(100, 90, 0.0, 0.03, 0.04, 1.5, 0.04, 0.5, -0.7,
                       1.0, -0.1, 0.2, OptionType.CALL) == pytest.approx(10.0)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        bates_price(100, 100, 1.0, 0.03, -0.04, 1.5, 0.04, 0.5, -0.7,
                    1.0, -0.1, 0.2)
    with pytest.raises(ValueError):
        bates_price(100, 100, 1.0, 0.03, 0.04, 1.5, 0.04, 0.5, -0.7,
                    -1.0, -0.1, 0.2)  # negative intensity


def test_deep_itm_call_approx_forward_minus_strike():
    # A deep in-the-money call is worth ~ discounted (forward - strike). The
    # fixed 64-node Gauss-Legendre integral leaves a small residual this deep, so
    # allow a slightly looser band than at-the-money.
    S, K, t, r, q = 100, 40, 1.0, 0.03, 0.0
    c = bates_price(S, K, t, r, 0.04, 1.5, 0.04, 0.5, -0.7,
                    1.0, -0.1, 0.2, OptionType.CALL, q=q)
    fwd = S * math.exp((r - q) * t)
    assert c == pytest.approx(math.exp(-r * t) * (fwd - K), abs=0.15)


@pytest.mark.slow
@pytest.mark.parametrize("S,K,t,r,q,args", [
    (100, 100, 0.5, 0.03, 0.0, (0.04, 1.5, 0.04, 0.5, -0.7, 1.0, -0.1, 0.15)),
    (100, 90, 1.0, 0.05, 0.02, (0.05, 1.0, 0.05, 0.6, -0.5, 2.0, -0.15, 0.2)),
    (100, 110, 0.25, 0.03, 0.0, (0.04, 2.0, 0.04, 0.4, -0.6, 3.0, -0.05, 0.1)),
])
def test_char_function_matches_qe_jump_monte_carlo(S, K, t, r, q, args):
    # Independent cross-check: an Andersen-QE variance simulation with a
    # compound-Poisson jump overlay must reproduce the Fourier price.
    import math as _m
    import random
    from quantforge.heston_mc import _norm_ppf

    v0, kappa, theta, xi, rho, lam, mu_j, sigma_j = args
    ns, npths, seed = 80, 150_000, 7
    dt = t / ns
    ekt = _m.exp(-kappa * dt)
    g1 = g2 = 0.5
    K0 = -rho * kappa * theta * dt / xi
    K1 = g1 * dt * (kappa * rho / xi - 0.5) - rho / xi
    K2 = g2 * dt * (kappa * rho / xi - 0.5) + rho / xi
    K3 = g1 * dt * (1 - rho * rho)
    K4 = g2 * dt * (1 - rho * rho)
    kbar = _m.exp(mu_j + 0.5 * sigma_j * sigma_j) - 1.0
    K0d = K0 + (r - q - lam * kbar) * dt
    rng = random.Random(seed)
    disc = _m.exp(-r * t)

    def nv(v, uz):
        m = theta + (v - theta) * ekt
        s2 = (v * xi * xi * ekt / kappa) * (1 - ekt) \
            + (theta * xi * xi / (2 * kappa)) * (1 - ekt) ** 2
        if m <= 0:
            return 0.0
        psi = s2 / (m * m)
        if psi <= 1.5:
            inv = 2 / psi
            b2 = inv - 1 + _m.sqrt(inv) * _m.sqrt(max(inv - 1, 0))
            a = m / (1 + b2)
            return a * (_m.sqrt(b2) + _norm_ppf(uz)) ** 2
        p = (psi - 1) / (psi + 1)
        beta = (1 - p) / m
        return 0.0 if uz <= p else _m.log((1 - p) / (1 - uz)) / beta

    sm = []
    for _ in range(npths // 2):
        z = [rng.gauss(0, 1) for _ in range(ns)]
        u = [rng.random() for _ in range(ns)]
        jz = [rng.gauss(0, 1) for _ in range(ns)]
        ju = [rng.random() for _ in range(ns)]
        for sgn in (1, -1):
            x = _m.log(S)
            v = v0
            for i in range(ns):
                vn = nv(v, u[i] if sgn == 1 else 1 - u[i])
                vol = _m.sqrt(max(K3 * v + K4 * vn, 0))
                x += K0d + K1 * v + K2 * vn + vol * (sgn * z[i])
                cnt, term, cum, uu = 0, _m.exp(-lam * dt), _m.exp(-lam * dt), ju[i]
                while uu > cum and cnt < 50:
                    cnt += 1
                    term *= lam * dt / cnt
                    cum += term
                if cnt > 0:
                    x += cnt * mu_j + _m.sqrt(cnt) * sigma_j * jz[i]
                v = vn
            sm.append(disc * max(_m.exp(x) - K, 0))
    mean = sum(sm) / len(sm)
    var = sum((a - mean) ** 2 for a in sm) / (len(sm) - 1)
    se = _m.sqrt(var / len(sm))

    cf = bates_price(S, K, t, r, *args, OptionType.CALL, q=q)
    assert abs(cf - mean) < 4.0 * se + 0.01
