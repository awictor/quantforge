"""Cohen's, weighted, and Fleiss' kappa."""

import random

import pytest

from quantforge import cohen_kappa, weighted_kappa, fleiss_kappa


def test_perfect_agreement():
    a = [1, 2, 3, 1, 2, 3]
    assert cohen_kappa(a, a) == 1.0


def test_cohen_hand_computed():
    a = ['y', 'y', 'n', 'n', 'y']
    b = ['y', 'n', 'n', 'n', 'y']
    # p_obs=0.8, p_exp=0.48 -> (0.8-0.48)/0.52 = 0.6154.
    assert abs(cohen_kappa(a, b) - 0.615385) < 1e-5


def test_chance_level_near_zero():
    rng = random.Random(3)
    a = [rng.randint(0, 2) for _ in range(3000)]
    b = [rng.randint(0, 2) for _ in range(3000)]
    assert abs(cohen_kappa(a, b)) < 0.05


def test_two_category_linear_weighted_equals_cohen():
    a = [0, 0, 1, 1, 0]
    b = [0, 1, 1, 1, 0]
    assert abs(cohen_kappa(a, b) - weighted_kappa(a, b, "linear")) < 1e-9


def test_weighted_rewards_near_misses():
    near = ([1, 2, 3, 2, 1], [1, 2, 2, 2, 1])   # off by 1
    far = ([1, 2, 3, 2, 1], [3, 2, 1, 2, 3])    # off by 2
    assert weighted_kappa(*near, "linear") > weighted_kappa(*far, "linear")


def test_fleiss_textbook():
    table = [[0, 0, 0, 0, 14], [0, 2, 6, 4, 2], [0, 0, 3, 5, 6],
             [0, 3, 9, 2, 0], [2, 2, 8, 1, 1]]
    assert abs(fleiss_kappa(table) - 0.2519) < 1e-3


def test_fleiss_perfect():
    assert fleiss_kappa([[3, 0], [0, 3], [3, 0]]) == 1.0


def test_validation():
    with pytest.raises(ValueError):
        cohen_kappa([1, 2], [1])
    with pytest.raises(ValueError):
        weighted_kappa([1, 2], [1, 2], weights="bad")
    with pytest.raises(ValueError):
        fleiss_kappa([[3, 0]])                        # one subject
    with pytest.raises(ValueError):
        fleiss_kappa([[3, 0], [2, 0]])                # unequal ratings per subject
