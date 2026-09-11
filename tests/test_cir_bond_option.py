"""European option on a CIR zero-coupon bond (cir_bond_option, CIR 1985)."""

import math
import random

import pytest

from quantforge import cir_bond_option, cir_zero_coupon_bond, OptionType
from quantforge.cir import _AB


KAPPA, THETA, SIG, R0 = 0.3, 0.04, 0.05, 0.04
S, T, K = 1.0, 3.0, 0.90


def test_put_call_parity():
    c = cir_bond_option(R0, S, T, K, KAPPA, THETA, SIG, OptionType.CALL)
    p = cir_bond_option(R0, S, T, K, KAPPA, THETA, SIG, OptionType.PUT)
    PS = cir_zero_coupon_bond(R0, S, KAPPA, THETA, SIG)
    PT = cir_zero_coupon_bond(R0, T, KAPPA, THETA, SIG)
    assert (c - p) == pytest.approx(PT - K * PS, abs=1e-9)


def test_prices_positive():
    c = cir_bond_option(R0, S, T, K, KAPPA, THETA, SIG, OptionType.CALL)
    p = cir_bond_option(R0, S, T, K, KAPPA, THETA, SIG, OptionType.PUT)
    assert c > 0.0 and p > 0.0


@pytest.mark.slow
def test_call_matches_monte_carlo():
    c = cir_bond_option(R0, S, T, K, KAPPA, THETA, SIG, OptionType.CALL)
    rng = random.Random(1)
    n = 200_000
    dt = 1.0 / 400
    nstep = int(S / dt)
    A_TS, B_TS = _AB(T - S, KAPPA, THETA, SIG)
    tot = 0.0
    for _ in range(n):
        r = R0
        integ = 0.0
        for _ in range(nstep):
            integ += r * dt
            r = max(r + KAPPA * (THETA - r) * dt
                    + SIG * math.sqrt(max(r, 0.0) * dt) * rng.gauss(0.0, 1.0), 0.0)
        PST = A_TS * math.exp(-B_TS * r)
        tot += math.exp(-integ) * max(PST - K, 0.0)
    assert c == pytest.approx(tot / n, abs=5e-4)


def test_higher_strike_lowers_call():
    lo = cir_bond_option(R0, S, T, 0.88, KAPPA, THETA, SIG, OptionType.CALL)
    hi = cir_bond_option(R0, S, T, 0.92, KAPPA, THETA, SIG, OptionType.CALL)
    assert hi < lo


def test_bad_times_raise():
    with pytest.raises(ValueError):
        cir_bond_option(R0, 3.0, 1.0, K, KAPPA, THETA, SIG)
