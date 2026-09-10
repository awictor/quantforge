"""Variance-swap term structure and forward-variance curve."""

import pytest

from quantforge import (
    variance_term_structure as vts,
    forward_variance_swap_from_smile as fwd,
)


S0, R = 100.0, 0.03
EXPIRIES = [0.25, 0.5, 1.0, 2.0]


def test_flat_term_structure_is_constant():
    flat = [lambda K: 0.2] * 4
    spot, forward = vts(S0, R, EXPIRIES, flat, n_strikes=101)
    # The strip's log-moneyness width scales with sqrt(t), so per-expiry
    # truncation differs at the ~1e-6 level; the curve is constant to that.
    assert all(s == pytest.approx(spot[0], abs=1e-4) for s in spot)
    assert all(f == pytest.approx(spot[0], abs=1e-4) for f in forward)


def test_rising_vol_gives_increasing_forward_variance():
    rising = [lambda K, s=s: s for s in (0.15, 0.18, 0.22, 0.26)]
    _spot, forward = vts(S0, R, EXPIRIES, rising, n_strikes=101)
    assert all(forward[i] <= forward[i + 1] + 1e-9 for i in range(3))
    assert all(f > 0 for f in forward)


def test_forward_curve_matches_forward_swap():
    rising = [lambda K, s=s: s for s in (0.15, 0.18, 0.22, 0.26)]
    _spot, forward = vts(S0, R, EXPIRIES, rising, n_strikes=101)
    fw = fwd(S0, EXPIRIES[1], EXPIRIES[2], R, rising[1], rising[2],
             n_strikes=101)
    assert forward[2] == pytest.approx(fw, abs=1e-9)


def test_first_forward_equals_first_spot():
    flat = [lambda K: 0.2] * 4
    spot, forward = vts(S0, R, EXPIRIES, flat, n_strikes=101)
    assert forward[0] == pytest.approx(spot[0])


def test_bad_inputs_raise():
    with pytest.raises(ValueError):
        vts(S0, R, [0.5, 1.0], [lambda K: 0.2])          # length mismatch
    with pytest.raises(ValueError):
        vts(S0, R, [1.0, 0.5], [lambda K: 0.2, lambda K: 0.2])  # not increasing
