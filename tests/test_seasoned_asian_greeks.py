"""Greeks of seasoned (in-progress) discrete Asian options (FD)."""

import pytest

from quantforge import (
    seasoned_geometric_asian, seasoned_geometric_asian_greeks,
    seasoned_arithmetic_asian, seasoned_arithmetic_asian_greeks,
    discrete_geometric_asian_greeks, discrete_arithmetic_asian_greeks,
)


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2
FIELDS = ("price", "delta", "gamma", "vega", "theta")


def test_geo_no_observations_matches_fresh_greeks():
    n = 12
    seas = seasoned_geometric_asian_greeks(S, K, T, R, SIG, [], n, option_type="call")
    fresh = discrete_geometric_asian_greeks(S, K, T, R, SIG, n_fixings=n,
                                            option_type="call")
    for f in FIELDS:
        assert seas[f] == pytest.approx(fresh[f], abs=1e-6)


def test_arith_no_observations_matches_fresh_greeks():
    n = 12
    seas = seasoned_arithmetic_asian_greeks(S, K, T, R, SIG, [], n, option_type="call")
    fresh = discrete_arithmetic_asian_greeks(S, K, T, R, SIG, n_fixings=n,
                                             option_type="call")
    for f in FIELDS:
        assert seas[f] == pytest.approx(fresh[f], abs=1e-6)


def test_geo_price_field_matches_pricer():
    t, n = 0.5, 12
    rem = [t * i / 6 for i in range(1, 7)]
    obs = [100.0] * 6
    g = seasoned_geometric_asian_greeks(S, K, t, R, SIG, obs, n, rem, "call")
    assert g["price"] == pytest.approx(
        seasoned_geometric_asian(S, K, t, R, SIG, obs, n, rem, "call"), abs=1e-9)


def test_arith_price_field_matches_pricer():
    t, n = 0.5, 12
    rem = [t * i / 6 for i in range(1, 7)]
    obs = [100.0] * 6
    g = seasoned_arithmetic_asian_greeks(S, K, t, R, SIG, obs, n, rem, "call")
    assert g["price"] == pytest.approx(
        seasoned_arithmetic_asian(S, K, t, R, SIG, obs, n, rem, "call"), abs=1e-9)


@pytest.mark.parametrize("fn", [seasoned_geometric_asian_greeks,
                                seasoned_arithmetic_asian_greeks])
def test_partial_greek_signs(fn):
    t, n = 0.5, 12
    rem = [t * i / 6 for i in range(1, 7)]
    g = fn(S, K, t, R, SIG, [100.0] * 6, n, rem, "call")
    assert g["delta"] > 0.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0
    assert g["theta"] < 0.0  # long call loses value as time passes here


def test_seasoning_reduces_delta_and_vega():
    # As more fixings lock in, less of the average is still stochastic, so both
    # delta (sensitivity to spot) and vega shrink versus a fresh option.
    t, n = 0.5, 12
    rem6 = [t * i / 6 for i in range(1, 7)]
    seasoned = seasoned_arithmetic_asian_greeks(S, K, t, R, SIG, [100.0] * 6, n,
                                                rem6, "call")
    remall = [t * i / n for i in range(1, n + 1)]
    fresh = seasoned_arithmetic_asian_greeks(S, K, t, R, SIG, [], n, remall, "call")
    assert seasoned["delta"] < fresh["delta"]
    assert seasoned["vega"] < fresh["vega"]
