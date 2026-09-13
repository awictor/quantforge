"""Mutual information and transfer entropy."""

import random

import pytest

from quantforge import mutual_information, transfer_entropy


def test_mutual_information_independent_near_zero():
    rng = random.Random(1)
    x = [rng.gauss(0, 1) for _ in range(5000)]
    y = [rng.gauss(0, 1) for _ in range(5000)]
    assert mutual_information(x, y) < 0.05


def test_mutual_information_symmetric():
    rng = random.Random(1)
    x = [rng.gauss(0, 1) for _ in range(3000)]
    y = [xi + 0.1 * rng.gauss(0, 1) for xi in x]
    assert abs(mutual_information(x, y) - mutual_information(y, x)) < 1e-12


def test_mutual_information_high_for_dependence():
    rng = random.Random(2)
    x = [rng.gauss(0, 1) for _ in range(4000)]
    y = [xi + 0.1 * rng.gauss(0, 1) for xi in x]
    assert mutual_information(x, y) > 0.5


def test_transfer_entropy_direction():
    rng = random.Random(2)
    x = [rng.gauss(0, 1) for _ in range(5000)]
    y = [0.0]
    for t in range(4999):
        y.append(0.6 * x[t] + 0.3 * rng.gauss(0, 1))    # Y_{t+1} driven by X_t
    te_xy = transfer_entropy(x, y)
    te_yx = transfer_entropy(y, x)
    assert te_xy > 3 * te_yx        # information flows X -> Y


def test_transfer_entropy_independent_near_zero():
    rng = random.Random(3)
    a = [rng.gauss(0, 1) for _ in range(5000)]
    b = [rng.gauss(0, 1) for _ in range(5000)]
    assert transfer_entropy(a, b) < 0.05


def test_non_negative():
    rng = random.Random(4)
    x = [rng.gauss(0, 1) for _ in range(2000)]
    y = [rng.gauss(0, 1) for _ in range(2000)]
    assert mutual_information(x, y) >= 0.0
    assert transfer_entropy(x, y) >= 0.0


def test_validation():
    with pytest.raises(ValueError):
        mutual_information([1.0], [1.0])
    with pytest.raises(ValueError):
        mutual_information([1.0, 2.0], [1.0, 2.0], bins=1)
    with pytest.raises(ValueError):
        transfer_entropy([1.0, 2.0], [1.0, 2.0])       # < 3 points
