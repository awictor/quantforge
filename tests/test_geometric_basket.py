"""Weighted geometric-average basket option (multiasset.geometric_basket_option)."""

import math

import pytest

from quantforge import geometric_basket_option, geometric_basket_greeks
from quantforge.bsm import price, OptionType


def test_single_asset_recovers_vanilla():
    for ot in (OptionType.CALL, OptionType.PUT):
        assert geometric_basket_option([100.0], [1.0], 100.0, 1.0, 0.05, [0.2],
                                       [[1.0]], None, ot) == pytest.approx(
            price(100.0, 100.0, 1.0, 0.05, 0.2, ot), abs=1e-10)


S = [100.0, 95.0, 105.0]
W = [0.5, 0.3, 0.2]
SIG = [0.2, 0.25, 0.18]
K, T, R = 100.0, 1.0, 0.05
CORR = [[1.0, 0.3, 0.2], [0.3, 1.0, 0.4], [0.2, 0.4, 1.0]]


def _fwd_var():
    mean = sum(W[i] * (math.log(S[i]) + (R - 0.5 * SIG[i] ** 2) * T) for i in range(3))
    var = T * sum(W[i] * W[j] * CORR[i][j] * SIG[i] * SIG[j]
                  for i in range(3) for j in range(3))
    return math.exp(mean + 0.5 * var)


def test_put_call_parity():
    c = geometric_basket_option(S, W, K, T, R, SIG, CORR, None, "call")
    p = geometric_basket_option(S, W, K, T, R, SIG, CORR, None, "put")
    F = _fwd_var()
    assert c - p == pytest.approx(math.exp(-R * T) * (F - K), abs=1e-9)


def test_matches_monte_carlo():
    import random
    cf = geometric_basket_option(S, W, K, T, R, SIG, CORR, None, "call")
    # Cholesky of the 3x3 correlation (lower-triangular), pure-Python.
    a = CORR
    L = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(i + 1):
            ssum = sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                L[i][j] = math.sqrt(a[i][i] - ssum)
            else:
                L[i][j] = (a[i][j] - ssum) / L[j][j]
    rng = random.Random(3)
    N = 200000
    disc = math.exp(-R * T)
    acc = 0.0
    for _ in range(N):
        z = [rng.gauss(0.0, 1.0) for _ in range(3)]
        corrz = [sum(L[i][k] * z[k] for k in range(3)) for i in range(3)]
        B = 1.0
        for i in range(3):
            ST = S[i] * math.exp((R - 0.5 * SIG[i] ** 2) * T
                                 + SIG[i] * math.sqrt(T) * corrz[i])
            B *= ST ** W[i]
        acc += max(B - K, 0.0)
    mc = disc * acc / N
    assert cf == pytest.approx(mc, rel=0.02)


def test_zero_correlation_and_positive():
    corr0 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    v = geometric_basket_option(S, W, K, T, R, SIG, corr0, None, "call")
    assert v > 0.0


def test_greeks_structure_and_signs():
    g = geometric_basket_greeks(S, W, K, T, R, SIG, CORR, None, "call")
    assert g["price"] == pytest.approx(
        geometric_basket_option(S, W, K, T, R, SIG, CORR, None, "call"), abs=1e-9)
    assert len(g["delta"]) == 3
    assert len(g["gamma"]) == 3
    assert all(d > 0.0 for d in g["delta"])  # call rises with each spot
    assert g["vega"] > 0.0


def test_validation():
    with pytest.raises(ValueError):
        geometric_basket_option(S, [0.5, 0.5], K, T, R, SIG, CORR)  # length mismatch
    with pytest.raises(ValueError):
        geometric_basket_option(S, W, -1.0, T, R, SIG, CORR)  # bad strike
    with pytest.raises(ValueError):
        geometric_basket_option(S, W, K, T, R, SIG, [[1.0, 0.3], [0.3, 1.0]])  # bad corr
