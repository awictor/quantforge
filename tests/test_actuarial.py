"""Life-contingent actuarial functions."""

import pytest

from quantforge import (
    survival_probabilities, life_annuity_due, term_insurance,
    whole_life_insurance, pure_endowment, endowment_insurance,
)


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


def test_validation():
    with pytest.raises(ValueError):
        survival_probabilities([1.5])
    with pytest.raises(ValueError):
        life_annuity_due(PX, -1.5)
