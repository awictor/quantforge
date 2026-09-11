"""Seasoned (in-progress) discrete arithmetic Asian (exotics.seasoned_arithmetic_asian)."""

import math

import pytest

from quantforge import seasoned_arithmetic_asian, discrete_arithmetic_asian
from quantforge.bsm import OptionType


S, K, R, SIG = 100.0, 100.0, 0.05, 0.2


def test_no_observations_matches_fresh():
    t, n = 1.0, 12
    fresh = discrete_arithmetic_asian(S, K, t, R, SIG, n_fixings=n, option_type="call")
    seas = seasoned_arithmetic_asian(S, K, t, R, SIG, [], n, option_type="call")
    assert seas == pytest.approx(fresh, abs=1e-10)


def test_no_observations_matches_fresh_put():
    t, n = 0.6, 10
    fresh = discrete_arithmetic_asian(S, K, t, R, SIG, n_fixings=n, option_type="put")
    seas = seasoned_arithmetic_asian(S, K, t, R, SIG, [], n, option_type="put")
    assert seas == pytest.approx(fresh, abs=1e-10)


def test_all_observed_is_deterministic():
    obs = [100.0, 101.0, 99.0, 102.0, 98.0, 100.0]
    A = sum(obs) / len(obs)
    t = 1e-6
    val = seasoned_arithmetic_asian(S, K, t, R, SIG, obs, len(obs), [], "call")
    assert val == pytest.approx(math.exp(-R * t) * max(A - K, 0.0), abs=1e-6)


def test_deep_in_the_money_is_expected_average_minus_strike():
    # Huge observed fixings force the shifted strike K' = nK - Q <= 0, so the
    # call is always in the money: value = disc*(E[A] - K), put worthless.
    t, n = 0.5, 12
    rem = [t * i / 6 for i in range(1, 7)]
    obs = [1000.0] * 6
    Q = sum(obs)
    EA = (Q + S * sum(math.exp(R * tj) for tj in rem)) / n
    call = seasoned_arithmetic_asian(S, K, t, R, SIG, obs, n, rem, "call")
    put = seasoned_arithmetic_asian(S, K, t, R, SIG, obs, n, rem, "put")
    assert call == pytest.approx(math.exp(-R * t) * (EA - K), abs=1e-6)
    assert put == pytest.approx(0.0, abs=1e-12)


def test_partial_matches_monte_carlo():
    t = 0.5
    obs = [100.0] * 6
    rem = [t * i / 6 for i in range(1, 7)]
    cf = seasoned_arithmetic_asian(S, K, t, R, SIG, obs, 12, rem, "call")
    import random
    rng = random.Random(9)
    N = 300000
    disc = math.exp(-R * t)
    Q = sum(obs)
    acc = 0.0
    for _ in range(N):
        logS = math.log(S)
        prevt = 0.0
        fs = 0.0
        for i in range(1, 7):
            ti = rem[i - 1]
            dt = ti - prevt
            prevt = ti
            z = rng.gauss(0.0, 1.0)
            logS += (R - 0.5 * SIG * SIG) * dt + SIG * math.sqrt(dt) * z
            fs += math.exp(logS)
        A = (Q + fs) / 12
        acc += max(A - K, 0.0)
    mc = disc * acc / N
    assert cf == pytest.approx(mc, rel=0.01)


def test_seasoned_arithmetic_above_geometric():
    from quantforge import seasoned_geometric_asian
    t, n = 0.5, 12
    rem = [t * i / 6 for i in range(1, 7)]
    obs = [100.0] * 6
    ar = seasoned_arithmetic_asian(S, K, t, R, SIG, obs, n, rem, "call")
    ge = seasoned_geometric_asian(S, K, t, R, SIG, obs, n, rem, "call")
    assert ar > ge


def test_validation():
    with pytest.raises(ValueError):
        seasoned_arithmetic_asian(S, K, 0.5, R, SIG, [100.0], 12,
                                  remaining_fixing_times=[0.1])  # wrong length
    with pytest.raises(ValueError):
        seasoned_arithmetic_asian(S, K, 0.5, R, SIG, [0.0], 12)  # bad price
    with pytest.raises(ValueError):
        seasoned_arithmetic_asian(S, K, 0.5, R, SIG, [100.0] * 13, 12)  # too many
