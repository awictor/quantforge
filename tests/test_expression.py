"""Tests for the shunting-yard expression evaluator, cross-checked against Python eval."""

import math
import random

import pytest

from quantforge.expression import tokenize, shunting_yard, eval_rpn, eval_expression


def _close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def _gen(depth, rng):
    if depth <= 0 or rng.random() < 0.35:
        return str(rng.randint(1, 9))
    op = rng.choice(["+", "-", "*", "/"])  # no ^: Python int** towers hang the reference
    e = f"({_gen(depth - 1, rng)}{op}{_gen(depth - 1, rng)})"
    if rng.random() < 0.2:
        e = "-" + e
    return e


def test_fuzz_vs_python_eval():
    rng = random.Random(231)
    tested = 0
    for _ in range(8000):
        expr = _gen(rng.randint(0, 4), rng)
        try:
            ref = eval(expr, {"__builtins__": {}}, {})
        except ZeroDivisionError:
            continue
        assert _close(eval_expression(expr), ref)
        tested += 1
    assert tested > 5000


def test_precedence():
    assert _close(eval_expression("2+3*4"), 14)
    assert _close(eval_expression("(2+3)*4"), 20)
    assert _close(eval_expression("2+3*4-1"), 13)
    assert _close(eval_expression("10-2-3"), 5)  # left-associative


def test_power_right_associative():
    assert _close(eval_expression("2^3^2"), 512)  # 2^(3^2) = 2^9
    assert _close(eval_expression("2^2^3"), 256)  # 2^(2^3) = 2^8
    assert _close(eval_expression("2^-1"), 0.5)


def test_unary_minus():
    assert _close(eval_expression("-2^2"), -4)   # -(2^2), power binds tighter
    assert _close(eval_expression("-(3+4)"), -7)
    assert _close(eval_expression("2*-3"), -6)
    assert _close(eval_expression("--3"), 3)


def test_division_variants():
    assert _close(eval_expression("10/4"), 2.5)
    assert _close(eval_expression("10//4"), 2)
    assert _close(eval_expression("10%3"), 1)


def test_functions_and_constants():
    assert _close(eval_expression("sqrt(16)"), 4)
    assert _close(eval_expression("sin(0)"), 0)
    assert _close(eval_expression("cos(pi)"), -1)
    assert _close(eval_expression("exp(1)"), math.e)
    assert _close(eval_expression("log(e)"), 1)
    assert _close(eval_expression("abs(-5)"), 5)
    assert _close(eval_expression("2*pi"), 2 * math.pi)


def test_nested_functions():
    assert _close(eval_expression("sqrt(3^2+4^2)"), 5)
    assert _close(eval_expression("abs(sin(0)-1)"), 1)


def test_scientific_notation():
    assert _close(eval_expression("1.5e3"), 1500)
    assert _close(eval_expression("2e-2"), 0.02)


def test_tokenize():
    assert tokenize("2 + 3") == [2.0, "+", 3.0]
    assert tokenize("(1+2)*3") == ["(", 1.0, "+", 2.0, ")", "*", 3.0]


def test_shunting_yard_and_eval_rpn_directly():
    rpn = shunting_yard(tokenize("2 + 3 * 4"))
    assert eval_rpn(rpn) == 14


def test_whitespace_insensitive():
    assert _close(eval_expression("  2  +  3  "), 5)


@pytest.mark.parametrize("bad", ["(1+2", "1+2)", "1++", "", "*3", "1 2", "@"])
def test_malformed_raises(bad):
    with pytest.raises((ValueError, IndexError)):
        eval_expression(bad)


def test_division_by_zero_propagates():
    with pytest.raises(ZeroDivisionError):
        eval_expression("1/0")
