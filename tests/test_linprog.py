"""Linear programming by two-phase simplex."""

import pytest

from quantforge import linprog


def test_max_with_le_constraints():
    r = linprog([3, 2], [([1, 1], "<=", 4), ([1, 3], "<=", 6)])
    assert abs(r["objective"] - 12.0) < 1e-6      # optimum at (4, 0)
    assert abs(r["x"][0] - 4.0) < 1e-6 and abs(r["x"][1]) < 1e-6


def test_box_constraints():
    r = linprog([1, 1], [([1, 0], "<=", 4), ([0, 1], "<=", 3)])
    assert abs(r["objective"] - 7.0) < 1e-6


def test_ge_constraint_phase1():
    r = linprog([2, 3], [([1, 1], ">=", 10)], maximize=False)
    assert abs(r["objective"] - 20.0) < 1e-6      # (10, 0)


def test_equality_constraint():
    r = linprog([1, 1], [([1, 1], "=", 5), ([1, 0], "<=", 3)])
    assert abs(r["objective"] - 5.0) < 1e-6


def test_diet_style_min_satisfies_constraints():
    r = linprog([0.6, 0.35],
                [([5, 7], ">=", 8), ([4, 2], ">=", 15), ([2, 1], ">=", 3)],
                maximize=False)
    x = r["x"]
    assert 5 * x[0] + 7 * x[1] >= 8 - 1e-6
    assert 4 * x[0] + 2 * x[1] >= 15 - 1e-6
    assert 2 * x[0] + x[1] >= 3 - 1e-6


def test_unbounded_and_infeasible():
    with pytest.raises(ValueError):
        linprog([1], [([1], ">=", 1)])                 # unbounded above
    with pytest.raises(ValueError):
        linprog([1], [([1], "<=", 1), ([1], ">=", 5)])  # infeasible
