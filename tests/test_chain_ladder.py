"""Chain-ladder loss reserving."""

import pytest

from quantforge import development_factors, chain_ladder


TRI = [[100, 150, 180], [110, 165], [120]]


def test_development_factors():
    f = development_factors(TRI)
    assert abs(f[0] - 1.5) < 1e-9      # (150+165)/(100+110)
    assert abs(f[1] - 1.2) < 1e-9      # 180/150


def test_ultimates_and_reserves():
    r = chain_ladder(TRI)
    assert r["ultimate"] == [180, 198.0, 216.0]
    assert r["reserve"] == [0, 33.0, 96.0]
    assert r["total_reserve"] == 129.0


def test_fully_developed_row_has_zero_reserve():
    assert chain_ladder(TRI)["reserve"][0] == 0


def test_factors_at_least_one_and_ultimate_above_latest():
    r = chain_ladder(TRI)
    assert all(x >= 1 for x in r["factors"])
    assert all(r["ultimate"][i] >= TRI[i][-1] for i in range(len(TRI)))


def test_validation():
    with pytest.raises(ValueError):
        development_factors([[100]])
