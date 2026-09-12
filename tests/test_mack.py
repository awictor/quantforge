"""Mack (1993) chain-ladder standard error, checked against the Taylor-Ashe triangle."""

import pytest

from quantforge import mack_standard_error, chain_ladder


# The cumulative run-off triangle from Mack (1993), ASTIN Bulletin 23(2), for which
# he reports a total reserve of 18,680,856 and a standard error of 2,447,095.
TAYLOR_ASHE = [
    [357848, 1124788, 1735330, 2218270, 2745596, 3319994, 3466336, 3606286, 3833515, 3901463],
    [352118, 1236139, 2170033, 3353322, 3799067, 4120063, 4647867, 4914039, 5339085],
    [290507, 1292306, 2218525, 3235179, 3985995, 4132918, 4628910, 4909315],
    [310608, 1418858, 2195047, 3757447, 4029929, 4381982, 4588268],
    [443160, 1136350, 2128333, 2897821, 3402672, 3873311],
    [396132, 1333217, 2180715, 2985752, 3691712],
    [440832, 1288463, 2419861, 3483130],
    [359480, 1421128, 2864498],
    [376686, 1363294],
    [344014],
]


def test_total_matches_mack_1993():
    r = mack_standard_error(TAYLOR_ASHE)
    assert round(r["total_reserve"]) == 18680856
    assert round(r["total_std_error"]) == 2447095


def test_per_year_std_errors_match_paper():
    r = mack_standard_error(TAYLOR_ASHE)
    # Mack's published per-accident-year standard errors (year 0 fully developed).
    expected = [0, 75535, 121699, 133549, 261406, 411010, 558317, 875328, 971258, 1363155]
    for i, e in enumerate(expected):
        assert abs(r["std_error"][i] - e) <= 1, (i, r["std_error"][i], e)


def test_ultimate_agrees_with_chain_ladder():
    r = mack_standard_error(TAYLOR_ASHE)
    cl = chain_ladder(TAYLOR_ASHE)
    for a, b in zip(r["ultimate"], cl["ultimate"]):
        assert abs(a - b) < 1e-6


def test_total_se_below_sum_of_year_se():
    # Positive between-year correlation makes the total SE exceed the root of the
    # summed variances but stay below the plain sum of the per-year SEs.
    r = mack_standard_error(TAYLOR_ASHE)
    quad = sum(s * s for s in r["std_error"]) ** 0.5
    linear = sum(r["std_error"])
    assert quad < r["total_std_error"] < linear


def test_cv_is_se_over_reserve():
    r = mack_standard_error(TAYLOR_ASHE)
    for i in range(len(TAYLOR_ASHE)):
        if r["reserve"][i] > 0:
            assert abs(r["cv"][i] - r["std_error"][i] / r["reserve"][i]) < 1e-12
    assert abs(r["total_cv"] - r["total_std_error"] / r["total_reserve"]) < 1e-12


def test_scale_invariance_of_cv():
    # Scaling every claim by a constant scales reserve and SE alike; the CV is fixed.
    scaled = [[2.0 * c for c in row] for row in TAYLOR_ASHE]
    a = mack_standard_error(TAYLOR_ASHE)
    b = mack_standard_error(scaled)
    assert abs(a["total_cv"] - b["total_cv"]) < 1e-10
    assert abs(b["total_reserve"] - 2.0 * a["total_reserve"]) < 1e-3


def test_validation():
    with pytest.raises(ValueError):
        mack_standard_error([[100.0]])
    with pytest.raises(ValueError):
        mack_standard_error([[100.0, -5.0], [200.0]])
