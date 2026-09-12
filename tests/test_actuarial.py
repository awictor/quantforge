"""Life-contingent actuarial functions."""

import pytest

import math

from quantforge import (
    survival_probabilities, life_annuity_due, term_insurance,
    whole_life_insurance, pure_endowment, endowment_insurance,
    temporary_life_annuity_due, net_level_premium,
    gompertz_makeham_hazard, gompertz_makeham_survival,
    gompertz_makeham_survival_curve, curtate_life_expectancy,
)


GM_A, GM_B, GM_C = 0.0005, 0.0000758, 1.09


PX = [0.99, 0.985, 0.98, 0.975, 0.97, 0.96, 0.95, 0.94, 0.92, 0.90]
I = 0.05


def test_survival_starts_at_one_and_decreases():
    cum = survival_probabilities(PX)
    assert cum[0] == 1.0
    assert all(cum[k] >= cum[k + 1] for k in range(len(cum) - 1))


def test_whole_life_annuity_identity():
    # A_x = 1 - d * a-due when the table terminates (everyone dead).
    px = PX + [0.0]
    adue = life_annuity_due(px, I)
    A = whole_life_insurance(px, I)
    d = I / (1 + I)
    assert A == pytest.approx(1 - d * adue, abs=1e-9)


def test_epvs_positive():
    assert life_annuity_due(PX, I) > 0
    assert term_insurance(PX, I) > 0


def test_term_below_whole_life():
    assert term_insurance(PX, I, 5) <= whole_life_insurance(PX, I) + 1e-12


def test_endowment_is_term_plus_pure_endowment():
    assert endowment_insurance(PX, I, 5) == pytest.approx(
        term_insurance(PX, I, 5) + pure_endowment(PX, I, 5), abs=1e-12)


def test_higher_interest_lowers_annuity():
    assert life_annuity_due(PX, 0.08) < life_annuity_due(PX, 0.05)


def test_pure_endowment_below_one():
    assert pure_endowment(PX, I, 5) < 1


def test_gm_hazard_increases_with_age():
    assert gompertz_makeham_hazard(70, GM_A, GM_B, GM_C) > \
        gompertz_makeham_hazard(40, GM_A, GM_B, GM_C)


def test_gm_survival_monotone_and_starts_at_one():
    svs = [gompertz_makeham_survival(40, t, GM_A, GM_B, GM_C) for t in range(20)]
    assert svs[0] == 1.0
    assert all(svs[i] >= svs[i + 1] for i in range(len(svs) - 1))


def test_gm_survival_matches_hazard_integral():
    def num_surv(x, t, n=10000):
        dt = t / n
        s = sum(gompertz_makeham_hazard(x + k * dt, GM_A, GM_B, GM_C) * dt
                for k in range(n))
        return math.exp(-s)
    assert gompertz_makeham_survival(40, 10, GM_A, GM_B, GM_C) == pytest.approx(
        num_surv(40, 10), abs=1e-4)


def test_gm_makeham_only_limit():
    # c -> 1: pure exponential (constant hazard a + b).
    assert gompertz_makeham_survival(40, 5, 0.01, 0.02, 1.0) == pytest.approx(
        math.exp(-(0.01 + 0.02) * 5), abs=1e-9)


def test_curtate_life_expectancy():
    px = gompertz_makeham_survival_curve(40, 60, GM_A, GM_B, GM_C)
    e = curtate_life_expectancy(px)
    assert e == pytest.approx(sum(survival_probabilities(px)[1:]), abs=1e-12)
    assert e > 0


def test_life_expectancy_decreases_with_age():
    e40 = curtate_life_expectancy(gompertz_makeham_survival_curve(40, 60, GM_A, GM_B, GM_C))
    e70 = curtate_life_expectancy(gompertz_makeham_survival_curve(70, 40, GM_A, GM_B, GM_C))
    assert e70 < e40


def test_gm_validation():
    with pytest.raises(ValueError):
        gompertz_makeham_survival(40, 5, GM_A, GM_B, 0)


def test_temporary_annuity_below_whole_life():
    assert temporary_life_annuity_due(PX, I, 5) <= life_annuity_due(PX, I) + 1e-12


def test_temporary_annuity_rises_with_term():
    assert temporary_life_annuity_due(PX, I, 3) < temporary_life_annuity_due(PX, I, 8)


def test_temporary_annuity_full_term_is_whole_life():
    assert temporary_life_annuity_due(PX, I, len(PX) + 1) == pytest.approx(
        life_annuity_due(PX, I), abs=1e-12)


def test_net_premium_equivalence_whole_life():
    P = net_level_premium(PX, I)
    assert P * life_annuity_due(PX, I) == pytest.approx(whole_life_insurance(PX, I),
                                                        abs=1e-12)
    assert 0 < P < 1


def test_net_premium_equivalence_endowment():
    P = net_level_premium(PX, I, 5)
    assert P * temporary_life_annuity_due(PX, I, 5) == pytest.approx(
        endowment_insurance(PX, I, 5), abs=1e-12)


def test_premium_validation():
    with pytest.raises(ValueError):
        temporary_life_annuity_due(PX, -1.5, 5)


def test_validation():
    with pytest.raises(ValueError):
        survival_probabilities([1.5])
    with pytest.raises(ValueError):
        life_annuity_due(PX, -1.5)
