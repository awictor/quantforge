"""Discount-only rho (carry fixed) for Black-76-style options (bsm module)."""

import pytest

from quantforge import rho_discount, price


S, K, T, SIG = 100.0, 105.0, 1.0, 0.2


@pytest.mark.parametrize("ot", ["call", "put"])
@pytest.mark.parametrize("r", [0.03, 0.05])
def test_matches_finite_difference_at_fixed_carry(ot, r):
    b = 0.0   # Black-76: carry independent of r
    h = 1e-6
    fd = (price(S, K, T, r + h, SIG, ot, b=b)
          - price(S, K, T, r - h, SIG, ot, b=b)) / (2 * h)
    assert rho_discount(S, K, T, r, SIG, ot, b=b) == pytest.approx(fd, abs=1e-4)


def test_equals_minus_t_times_price():
    r, b = 0.05, 0.0
    for ot in ("call", "put"):
        assert rho_discount(S, K, T, r, SIG, ot, b=b) == pytest.approx(
            -T * price(S, K, T, r, SIG, ot, b=b), abs=1e-12)


def test_both_negative():
    r, b = 0.05, 0.0
    assert rho_discount(S, K, T, r, SIG, "call", b=b) < 0.0
    assert rho_discount(S, K, T, r, SIG, "put", b=b) < 0.0


def test_zero_time_is_zero():
    assert rho_discount(S, K, 0.0, r=0.05, sigma=SIG, option_type="call", b=0.0) == 0.0
