"""Multi-factor OLS return regression."""

import pytest

from quantforge import factor_regression, factor_expected_return
from quantforge.correlation import realized_beta


F1 = [0.02, -0.01, 0.03, 0.00, 0.015, -0.02, 0.025, 0.01]
F2 = [0.01, 0.02, -0.01, 0.015, 0.00, 0.01, -0.02, 0.005]
R = [0.01 + 1.5 * a - 0.3 * b for a, b in zip(F1, F2)]


def test_recovers_exact_coefficients():
    res = factor_regression(R, [F1, F2])
    assert res["alpha"] == pytest.approx(0.01, abs=1e-9)
    assert res["betas"][0] == pytest.approx(1.5, abs=1e-9)
    assert res["betas"][1] == pytest.approx(-0.3, abs=1e-9)


def test_exact_fit_r_squared_one():
    res = factor_regression(R, [F1, F2])
    assert res["r_squared"] == pytest.approx(1.0, abs=1e-9)
    assert res["residual_vol"] < 1e-9


def test_single_factor_matches_realized_beta():
    res = factor_regression(R, [F1])
    assert res["betas"][0] == pytest.approx(realized_beta(R, F1), abs=1e-9)


def test_expected_return():
    assert factor_expected_return(0.01, [1.5, -0.3], [0.06, 0.02]) == pytest.approx(
        0.01 + 1.5 * 0.06 - 0.3 * 0.02)


def test_noisy_fit_r_squared_below_one():
    rn = [R[i] + (0.001 if i % 2 else -0.001) for i in range(8)]
    assert factor_regression(rn, [F1, F2])["r_squared"] < 1.0


def test_validation():
    with pytest.raises(ValueError):
        factor_regression(R, [F1[:3]])
    with pytest.raises(ValueError):
        factor_expected_return(0.01, [1.5], [0.06, 0.02])
