"""GARCH term-vol option pricing bridge (volatility.garch_option_price)."""

import math

import pytest

from quantforge import (
    GarchParams, garch_option_price, garch_term_variance,
)
from quantforge.bsm import price


P = GarchParams(omega=2e-6, alpha=0.08, beta=0.90)
LR = P.long_run_variance
H = LR * 1.5
R = 0.01
S, K, RR = 100.0, 100.0, 0.05
N = 21   # ~1 trading month


@pytest.mark.parametrize("ot", ["call", "put"])
def test_matches_bsm_with_term_vol(ot):
    sig = garch_term_variance(P, R, H, N, 252)
    t = N / 252
    assert garch_option_price(P, R, H, S, K, RR, ot, horizon=N) == pytest.approx(
        price(S, K, t, RR, sig, ot), abs=1e-12)


def test_elevated_variance_richer_than_long_run():
    lr_sig = (LR * 252) ** 0.5
    g = garch_option_price(P, R, H, S, K, RR, "call", horizon=5)
    lr_price = price(S, K, 5 / 252, RR, lr_sig, "call")
    assert g > lr_price


def test_explicit_t_override():
    sig = garch_term_variance(P, R, H, N, 252)
    assert garch_option_price(P, R, H, S, K, RR, "call", horizon=N, t=0.5) == pytest.approx(
        price(S, K, 0.5, RR, sig, "call"), abs=1e-12)


def test_put_call_parity():
    c = garch_option_price(P, R, H, S, K, RR, "call", horizon=N)
    p = garch_option_price(P, R, H, S, K, RR, "put", horizon=N)
    t = N / 252
    assert c - p == pytest.approx(S - K * math.exp(-RR * t), abs=1e-9)


def test_missing_horizon_raises():
    with pytest.raises(ValueError):
        garch_option_price(P, R, H, S, K, RR, "call")
