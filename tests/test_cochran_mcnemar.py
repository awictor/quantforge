"""McNemar and Cochran's Q tests for paired binary data."""

import random

import pytest

from quantforge import mcnemar_test, cochran_q_test
from quantforge.distributions import binomial_cdf


def test_mcnemar_exact_binomial():
    r = mcnemar_test(b=10, c=2)
    assert abs(r["p_value"] - min(2 * binomial_cdf(2, 12, 0.5), 1.0)) < 1e-12
    assert abs(r["chi2_cc"] - (abs(10 - 2) - 1) ** 2 / 12) < 1e-12


def test_mcnemar_table_matches_counts():
    r = mcnemar_test(table=[[5, 10], [2, 8]])
    assert r["b"] == 10 and r["c"] == 2


def test_mcnemar_equal_discordant_not_significant():
    assert mcnemar_test(b=7, c=7)["p_value"] == 1.0


def test_mcnemar_zero_discordant():
    r = mcnemar_test(b=0, c=0)
    assert r["p_value"] == 1.0 and r["chi2_cc"] == 0.0


def test_cochran_q_equals_uncorrected_mcnemar_for_two():
    rng = random.Random(3)
    rows = [[rng.randint(0, 1), rng.randint(0, 1)] for _ in range(50)]
    q = cochran_q_test(rows)["statistic"]
    b = sum(1 for r in rows if r[0] == 0 and r[1] == 1)
    c = sum(1 for r in rows if r[0] == 1 and r[1] == 0)
    mc = (b - c) ** 2 / (b + c) if b + c else 0
    assert abs(q - mc) < 1e-9


def test_cochran_q_detects_difference():
    data = [[0, 0, 1], [0, 1, 1], [0, 0, 1], [1, 0, 1], [0, 0, 1], [0, 1, 1]]
    r = cochran_q_test(data)
    assert r["df"] == 2
    assert r["p_value"] < 0.05


def test_cochran_q_null_calibration():
    rng = random.Random(9)
    rej = 0
    N = 1000
    for _ in range(N):
        rows = [[rng.randint(0, 1) for _ in range(3)] for _ in range(15)]
        if cochran_q_test(rows)["p_value"] < 0.05:
            rej += 1
    assert 0.03 <= rej / N <= 0.07


def test_validation():
    with pytest.raises(ValueError):
        mcnemar_test(table=[[1, 2, 3], [4, 5, 6]])
    with pytest.raises(ValueError):
        cochran_q_test([[0, 1, 1]])                    # one block
    with pytest.raises(ValueError):
        cochran_q_test([[0, 1], [1, 2]])               # non-binary
