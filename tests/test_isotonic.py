"""Isotonic regression by pool-adjacent-violators."""

import random

import pytest

from quantforge import isotonic_regression, isotonic_fit


def _minmax(y, w):
    # Independent reference: yhat_k = max_{i<=k} min_{j>=k} weighted_mean(y[i..j]).
    n = len(y)
    out = []
    for k in range(n):
        best = -1e18
        for i in range(k + 1):
            worst = 1e18
            for j in range(k, n):
                num = sum(w[t] * y[t] for t in range(i, j + 1))
                den = sum(w[t] for t in range(i, j + 1))
                worst = min(worst, num / den)
            best = max(best, worst)
        out.append(best)
    return out


def test_monotone_unchanged():
    assert isotonic_regression([1, 2, 3, 4, 5]) == [1.0, 2.0, 3.0, 4.0, 5.0]


def test_textbook_pooling():
    assert isotonic_regression([1, 2, 4, 2, 5]) == [1.0, 2.0, 3.0, 3.0, 5.0]


def test_matches_minmax_formula_weighted():
    rng = random.Random(3)
    for _ in range(300):
        n = rng.randint(1, 12)
        y = [rng.uniform(-5, 5) for _ in range(n)]
        w = [rng.uniform(0.1, 3) for _ in range(n)]
        a = isotonic_regression(y, w)
        b = _minmax(y, w)
        assert max(abs(a[i] - b[i]) for i in range(n)) < 1e-9


def test_result_is_monotone_and_optimal():
    rng = random.Random(9)
    n = 10
    y = [rng.uniform(0, 10) for _ in range(n)]
    fit = isotonic_regression(y)
    assert all(fit[i] <= fit[i + 1] + 1e-12 for i in range(n - 1))
    sse = lambda f: sum((f[i] - y[i]) ** 2 for i in range(n))
    # Beats sorted(y), which is another monotone competitor.
    assert sse(fit) <= sse(sorted(y)) + 1e-9


def test_decreasing():
    assert isotonic_regression([5, 2, 4, 2, 1], increasing=False) == [5.0, 3.0, 3.0, 2.0, 1.0]


def test_predict_clamps_and_interpolates():
    y_hat, predict = isotonic_fit([1, 2, 3, 4, 5], [1, 2, 4, 2, 5])
    assert predict(0) == 1.0        # clamp below
    assert predict(99) == 5.0       # clamp above
    assert abs(predict(3) - 3.0) < 1e-9

def test_fit_realigns_to_original_order():
    x = [3, 1, 2]
    y = [5, 1, 9]
    y_hat, _ = isotonic_fit(x, y)
    # sorted by x: (1,1),(2,9),(3,5) -> isotonic [1,7,7] -> realign to x order [7,1,7]
    assert abs(y_hat[0] - 7.0) < 1e-9
    assert abs(y_hat[1] - 1.0) < 1e-9
    assert abs(y_hat[2] - 7.0) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        isotonic_regression([1, 2], weights=[1])
    with pytest.raises(ValueError):
        isotonic_regression([1, 2], weights=[1, -1])
    with pytest.raises(ValueError):
        isotonic_fit([1, 2], [1])
