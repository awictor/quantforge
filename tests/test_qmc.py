"""Tests for the quasi-Monte Carlo (Halton) sequence and QMC pricer."""

import math

import pytest

from quantforge import halton, european_qmc, call_price, put_price, OptionType
from quantforge.qmc import _radical_inverse


def test_radical_inverse_base2_known_values():
    # 1 -> 0.1b = 0.5, 2 -> 0.01b = 0.25, 3 -> 0.11b = 0.75.
    assert _radical_inverse(1, 2) == pytest.approx(0.5)
    assert _radical_inverse(2, 2) == pytest.approx(0.25)
    assert _radical_inverse(3, 2) == pytest.approx(0.75)


def test_halton_points_in_unit_cube():
    for i in range(50):
        for x in halton(i, 3):
            assert 0.0 <= x < 1.0


def test_halton_is_deterministic():
    assert halton(10, 2) == halton(10, 2)


def test_qmc_converges_to_bsm_call():
    exact = call_price(100, 100, 1.0, 0.05, 0.2)
    assert european_qmc(100, 100, 1.0, 0.05, 0.2, OptionType.CALL, n_points=8192) \
        == pytest.approx(exact, abs=1e-2)


def test_qmc_converges_to_bsm_put():
    exact = put_price(100, 105, 0.5, 0.03, 0.3)
    assert european_qmc(100, 105, 0.5, 0.03, 0.3, OptionType.PUT, n_points=8192) \
        == pytest.approx(exact, abs=1e-2)


def test_qmc_error_shrinks_with_points():
    exact = call_price(100, 100, 1.0, 0.05, 0.2)
    coarse = abs(european_qmc(100, 100, 1.0, 0.05, 0.2, n_points=256) - exact)
    fine = abs(european_qmc(100, 100, 1.0, 0.05, 0.2, n_points=8192) - exact)
    assert fine < coarse


def test_qmc_deterministic_no_seed():
    a = european_qmc(100, 100, 1.0, 0.05, 0.2, n_points=1000)
    b = european_qmc(100, 100, 1.0, 0.05, 0.2, n_points=1000)
    assert a == b


def test_qmc_beats_pseudo_mc_on_average():
    # Averaged over several strikes, QMC error should be below pseudo-MC error
    # at the same point count.
    from quantforge import european_mc
    n = 2048
    q_err = m_err = 0.0
    for K in (80, 90, 100, 110, 120):
        exact = call_price(100, K, 1.0, 0.05, 0.2)
        q_err += abs(european_qmc(100, K, 1.0, 0.05, 0.2, n_points=n) - exact)
        m_err += abs(european_mc(100, K, 1.0, 0.05, 0.2, n_paths=n, seed=1).price - exact)
    assert q_err < m_err


def test_zero_vol_is_intrinsic():
    v = european_qmc(110, 100, 1.0, 0.0, 0.0, OptionType.CALL)
    assert v == pytest.approx(10.0, abs=1e-9)
