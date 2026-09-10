"""Two-asset correlation option: a vanilla on S1 gated by S2 (correlation_option)."""

import math
import random

import pytest

from quantforge import (
    correlation_option,
    call_price,
    put_price,
    OptionType,
)


S1, S2, K1, K2, T, R = 100.0, 95.0, 100.0, 90.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


def _mc(ot, cond2, seed, n=800_000):
    rng = random.Random(seed)
    d1 = (R - 0.5 * SIG1 * SIG1) * T
    d2 = (R - 0.5 * SIG2 * SIG2) * T
    v1 = SIG1 * math.sqrt(T)
    v2 = SIG2 * math.sqrt(T)
    c2 = math.sqrt(1.0 - RHO * RHO)
    disc = math.exp(-R * T)
    tot = 0.0
    for _ in range(n):
        z1 = rng.gauss(0.0, 1.0)
        z2 = rng.gauss(0.0, 1.0)
        a = S1 * math.exp(d1 + v1 * z1)
        b = S2 * math.exp(d2 + v2 * (RHO * z1 + c2 * z2))
        gate = (b > K2) if cond2 == "above" else (b < K2)
        if gate:
            tot += (max(a - K1, 0.0) if ot is OptionType.CALL
                    else max(K1 - a, 0.0))
    return disc * tot / n


@pytest.mark.slow
def test_call_matches_mc():
    closed = correlation_option(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                OptionType.CALL, "above")
    assert closed == pytest.approx(_mc(OptionType.CALL, "above", 1), abs=0.02)


@pytest.mark.slow
def test_put_matches_mc():
    closed = correlation_option(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                OptionType.PUT, "above")
    assert closed == pytest.approx(_mc(OptionType.PUT, "above", 2), abs=0.02)


@pytest.mark.slow
def test_below_barrier_matches_mc():
    closed = correlation_option(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                OptionType.CALL, "below")
    assert closed == pytest.approx(_mc(OptionType.CALL, "below", 3), abs=0.02)


def test_unconditional_reduces_to_vanilla_call():
    # A barrier that is essentially always met -> plain vanilla on S1.
    closed = correlation_option(S1, S2, K1, 1e-6, T, R, SIG1, SIG2, RHO,
                                OptionType.CALL, "above")
    assert closed == pytest.approx(call_price(S1, K1, T, R, SIG1), abs=1e-6)


def test_unconditional_reduces_to_vanilla_put():
    closed = correlation_option(S1, S2, K1, 1e-6, T, R, SIG1, SIG2, RHO,
                                OptionType.PUT, "above")
    assert closed == pytest.approx(put_price(S1, K1, T, R, SIG1), abs=1e-6)


def test_above_plus_below_is_vanilla():
    # Gating on S2>K2 and S2<K2 partitions all paths -> the two correlation
    # calls sum to the plain vanilla on S1.
    up = correlation_option(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                            OptionType.CALL, "above")
    dn = correlation_option(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                            OptionType.CALL, "below")
    assert up + dn == pytest.approx(call_price(S1, K1, T, R, SIG1), abs=1e-9)


def test_bad_cond_raises():
    with pytest.raises(ValueError):
        correlation_option(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO, cond2="maybe")


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        correlation_option(S1, S2, K1, K2, T, R, SIG1, SIG2, 1.5)
