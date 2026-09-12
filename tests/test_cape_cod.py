"""Cape Cod (Stanard-Buhlmann) reserving."""

import pytest

from quantforge import cape_cod, bornhuetter_ferguson
from quantforge.chain_ladder import development_factors, development_pattern


TRI = [[100, 150, 180], [110, 165], [120]]
PREM = [300, 300, 300]


def _pct():
    pat = development_pattern(development_factors(TRI))
    return [pat[len(TRI[i]) - 1] for i in range(len(TRI))]


def test_elr_is_loss_over_used_premium():
    pct = _pct()
    used = sum(PREM[i] * pct[i] for i in range(len(TRI)))
    expected = sum(TRI[i][-1] for i in range(len(TRI))) / used
    assert abs(cape_cod(TRI, PREM)["elr"] - expected) < 1e-6


def test_developed_year_zero_reserve():
    assert abs(cape_cod(TRI, PREM)["reserve"][0]) < 1e-9


def test_reserve_is_premium_times_elr_times_undeveloped():
    r = cape_cod(TRI, PREM)
    pct = _pct()
    assert abs(r["reserve"][2] - 300 * r["elr"] * (1 - pct[2])) < 1e-6


def test_equals_bf_with_premium_times_elr_apriori():
    r = cape_cod(TRI, PREM)
    bf = bornhuetter_ferguson(TRI, [PREM[i] * r["elr"] for i in range(len(TRI))])
    assert all(abs(r["reserve"][i] - bf["reserve"][i]) < 1e-6 for i in range(len(TRI)))


def test_validation():
    with pytest.raises(ValueError):
        cape_cod(TRI, [300, 300])
