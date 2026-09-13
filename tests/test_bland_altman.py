"""Bland-Altman agreement and Lin's concordance correlation."""

import math
import random
import statistics

import pytest

from quantforge import bland_altman, concordance_correlation


def test_identical_methods():
    x = [1, 2, 3, 4, 5]
    r = bland_altman(x, x)
    assert r["bias"] == 0.0
    assert r["sd"] == 0.0
    assert abs(concordance_correlation(x, x) - 1.0) < 1e-12


def test_bias_and_sd_hand_computed():
    x = [10, 20, 30, 40]
    y = [12, 19, 33, 38]
    diffs = [x[i] - y[i] for i in range(4)]
    r = bland_altman(x, y)
    assert abs(r["bias"] - sum(diffs) / 4) < 1e-12
    assert abs(r["sd"] - statistics.stdev(diffs)) < 1e-12
    assert abs(r["lower"] - (r["bias"] - 1.96 * r["sd"])) < 1e-12
    assert abs(r["upper"] - (r["bias"] + 1.96 * r["sd"])) < 1e-12


def test_ccc_below_pearson_under_offset():
    x = [1, 2, 3, 4, 5]
    y = [xi + 2 for xi in x]           # perfect correlation, constant offset
    # Pearson is 1, but CCC is penalized for the bias.
    assert abs(concordance_correlation(x, y) - 0.5) < 1e-9


def test_ccc_matches_reference():
    def ref(x, y):
        n = len(x)
        mx, my = sum(x) / n, sum(y) / n
        sx2 = sum((v - mx) ** 2 for v in x) / n
        sy2 = sum((v - my) ** 2 for v in y) / n
        sxy = sum((x[i] - mx) * (y[i] - my) for i in range(n)) / n
        return 2 * sxy / (sx2 + sy2 + (mx - my) ** 2)

    rng = random.Random(3)
    for _ in range(300):
        n = rng.randint(2, 20)
        xs = [rng.gauss(0, 2) for _ in range(n)]
        ys = [rng.gauss(0.5, 2) for _ in range(n)]
        assert abs(concordance_correlation(xs, ys) - ref(xs, ys)) < 1e-9


def test_ccc_never_exceeds_pearson():
    rng = random.Random(7)
    for _ in range(300):
        n = rng.randint(3, 20)
        xs = [rng.gauss(0, 2) for _ in range(n)]
        ys = [1.3 * v + rng.gauss(0.5, 1) for v in xs]
        mx, my = sum(xs) / n, sum(ys) / n
        denom = math.sqrt(sum((v - mx) ** 2 for v in xs) * sum((v - my) ** 2 for v in ys))
        if denom == 0:
            continue
        pearson = sum((xs[i] - mx) * (ys[i] - my) for i in range(n)) / denom
        assert abs(concordance_correlation(xs, ys)) <= abs(pearson) + 1e-9


def test_limits_cover_about_95_percent():
    rng = random.Random(9)
    x = [rng.gauss(50, 10) for _ in range(10000)]
    y = [xi + rng.gauss(0, 2) for xi in x]
    r = bland_altman(x, y)
    inside = sum(1 for d in r["diffs"] if r["lower"] <= d <= r["upper"]) / len(r["diffs"])
    assert 0.93 <= inside <= 0.97


def test_validation():
    with pytest.raises(ValueError):
        bland_altman([1, 2], [1])
    with pytest.raises(ValueError):
        concordance_correlation([1], [1])
