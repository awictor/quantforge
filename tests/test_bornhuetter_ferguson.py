"""Bornhuetter-Ferguson reserving."""

import pytest

from quantforge import (
    bornhuetter_ferguson, development_pattern, development_factors, chain_ladder,
)


TRI = [[100, 150, 180], [110, 165], [120]]


def test_development_pattern():
    pat = development_pattern(development_factors(TRI))
    assert abs(pat[2] - 1.0) < 1e-9
    assert abs(pat[1] - 1 / 1.2) < 1e-6
    assert abs(pat[0] - 1 / 1.8) < 1e-6


def test_reserve_is_apriori_times_undeveloped():
    r = bornhuetter_ferguson(TRI, [200, 200, 200])
    assert abs(r["reserve"][0]) < 1e-9                       # fully developed
    assert abs(r["reserve"][1] - 200 * (1 - 1 / 1.2)) < 1e-6
    assert abs(r["reserve"][2] - 200 * (1 - 1 / 1.8)) < 1e-6


def test_ultimate_is_latest_plus_reserve():
    r = bornhuetter_ferguson(TRI, [200, 200, 200])
    assert abs(r["ultimate"][2] - (120 + r["reserve"][2])) < 1e-9


def test_bf_with_cl_apriori_matches_chain_ladder():
    cl = chain_ladder(TRI)
    bf = bornhuetter_ferguson(TRI, cl["ultimate"])
    assert all(abs(bf["reserve"][i] - cl["reserve"][i]) < 1e-6 for i in range(len(TRI)))


def test_validation():
    with pytest.raises(ValueError):
        bornhuetter_ferguson(TRI, [200, 200])
