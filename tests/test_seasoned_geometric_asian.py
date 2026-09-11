"""Seasoned (in-progress) discrete geometric Asian (exotics.seasoned_geometric_asian)."""

import math

import pytest

from quantforge import seasoned_geometric_asian, discrete_geometric_asian
from quantforge.bsm import OptionType


S, K, R, SIG = 100.0, 100.0, 0.05, 0.2


def test_no_observations_matches_fresh():
    t, n = 1.0, 12
    fresh = discrete_geometric_asian(S, K, t, R, SIG, n_fixings=n, option_type="call")
    seas = seasoned_geometric_asian(S, K, t, R, SIG, [], n, option_type="call")
    assert seas == pytest.approx(fresh, abs=1e-10)


def test_no_observations_matches_fresh_put():
    t, n = 0.7, 8
    fresh = discrete_geometric_asian(S, K, t, R, SIG, n_fixings=n, option_type="put")
    seas = seasoned_geometric_asian(S, K, t, R, SIG, [], n, option_type="put")
    assert seas == pytest.approx(fresh, abs=1e-10)


def test_all_observed_is_deterministic():
    obs = [100.0, 101.0, 99.0, 102.0, 98.0, 100.0]
    G = math.exp(sum(math.log(x) for x in obs) / len(obs))
    t = 1e-6
    disc = math.exp(-R * t)
    val = seasoned_geometric_asian(S, K, t, R, SIG, obs, len(obs), [], "call")
    assert val == pytest.approx(disc * max(G - K, 0.0), abs=1e-6)


def test_all_observed_itm_call():
    obs = [110.0, 112.0, 108.0, 111.0]
    G = math.exp(sum(math.log(x) for x in obs) / len(obs))
    t = 1e-6
    val = seasoned_geometric_asian(S, K, t, R, SIG, obs, len(obs), [], "call")
    assert val == pytest.approx(math.exp(-R * t) * (G - K), abs=1e-4)


def test_partial_matches_monte_carlo():
    t = 0.5
    obs = [100.0] * 6
    rem = [t * i / 6 for i in range(1, 7)]
    cf = seasoned_geometric_asian(S, K, t, R, SIG, obs, 12, rem, "call")
    import random
    rng = random.Random(7)
    N = 300000
    disc = math.exp(-R * t)
    Aobs = sum(math.log(x) for x in obs)
    acc = 0.0
    for _ in range(N):
        logS = math.log(S)
        prevt = 0.0
        fut = 0.0
        for i in range(1, 7):
            ti = rem[i - 1]
            dt = ti - prevt
            prevt = ti
            z = rng.gauss(0.0, 1.0)
            logS += (R - 0.5 * SIG * SIG) * dt + SIG * math.sqrt(dt) * z
            fut += logS
        G = math.exp((Aobs + fut) / 12)
        acc += max(G - K, 0.0)
    mc = disc * acc / N
    assert cf == pytest.approx(mc, rel=0.01)


def test_observed_above_strike_raises_floor():
    # High observed fixings lock in part of the average -> call worth more than
    # the same option with at-the-money observations.
    t, n = 0.5, 12
    rem = [t * i / 6 for i in range(1, 7)]
    high = seasoned_geometric_asian(S, K, t, R, SIG, [130.0] * 6, n, rem, "call")
    atm = seasoned_geometric_asian(S, K, t, R, SIG, [100.0] * 6, n, rem, "call")
    assert high > atm


def test_validation():
    with pytest.raises(ValueError):
        seasoned_geometric_asian(S, K, 0.5, R, SIG, [100.0], 12,
                                 remaining_fixing_times=[0.1, 0.2])  # wrong length
    with pytest.raises(ValueError):
        seasoned_geometric_asian(S, K, 0.5, R, SIG, [-1.0], 12)  # bad price
    with pytest.raises(ValueError):
        seasoned_geometric_asian(S, K, 0.5, R, SIG, [100.0] * 13, 12)  # too many
