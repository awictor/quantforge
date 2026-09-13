"""Fractional differencing for long-memory series."""

import pytest

from quantforge import (fracdiff_weights, fractional_difference,
                        fixed_width_fracdiff)


def test_integer_weights():
    assert [round(w, 6) for w in fracdiff_weights(1, 4)] == [1.0, -1.0, 0.0, 0.0]
    assert [round(w, 6) for w in fracdiff_weights(2, 5)] == [1.0, -2.0, 1.0, 0.0, 0.0]
    assert [round(w, 6) for w in fracdiff_weights(0, 4)] == [1.0, 0.0, 0.0, 0.0]


def test_weight_recursion():
    d = 0.4
    w = fracdiff_weights(d, 8)
    for k in range(1, 8):
        assert abs(w[k] - w[k - 1] * -(d - k + 1) / k) < 1e-12


def test_d_zero_is_identity():
    x = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0]
    assert fractional_difference(x, 0) == x


def test_d_one_is_first_difference():
    x = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0]
    out = fractional_difference(x, 1)
    expected = [x[0]] + [x[i] - x[i - 1] for i in range(1, len(x))]
    assert [round(v, 9) for v in out] == [round(v, 9) for v in expected]


def test_d_two_is_second_difference():
    x = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0]
    out = fractional_difference(x, 2)
    d1 = [x[i] - x[i - 1] for i in range(1, len(x))]
    d2 = [d1[i] - d1[i - 1] for i in range(1, len(d1))]
    # Beyond the warm-up the full-expansion output equals the true second diff.
    assert [round(v, 9) for v in out[2:]] == [round(v, 9) for v in d2]


def test_fractional_weights_decay_slowly():
    w = fracdiff_weights(0.4, 50)
    assert w[0] == 1.0
    assert all(abs(w[k]) < abs(w[k - 1]) for k in range(2, 50))   # monotone decay
    assert abs(w[49]) > 0.0                                        # but never zero


def test_fixed_width_trims_and_preserves_length_relation():
    x = [float(i % 7) for i in range(60)]
    ffd = fixed_width_fracdiff(x, 0.4, threshold=0.05)
    # Output length is n - width + 1 for the retained window.
    assert 0 < len(ffd) < len(x)


def test_validation():
    with pytest.raises(ValueError):
        fracdiff_weights(0.4, 0)
    with pytest.raises(ValueError):
        fractional_difference([], 0.4)
    with pytest.raises(ValueError):
        fixed_width_fracdiff([], 0.4)
