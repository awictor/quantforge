"""Third-order vega Greek ultima = d(vomma)/d(sigma) (greeks2.ultima)."""

import pytest

from quantforge import ultima, vomma
from quantforge.bsm import price, OptionType


CASES = [
    (100.0, 100.0, 1.0, 0.05, 0.2),
    (100.0, 90.0, 0.5, 0.03, 0.3),
    (100.0, 110.0, 2.0, 0.02, 0.25),
    (50.0, 55.0, 0.75, 0.04, 0.35),
]


@pytest.mark.parametrize("S,K,t,r,sig", CASES)
def test_matches_fd_of_vomma(S, K, t, r, sig):
    h = 1e-4
    fd = (vomma(S, K, t, r, sig + h) - vomma(S, K, t, r, sig - h)) / (2 * h)
    assert ultima(S, K, t, r, sig) == pytest.approx(fd, rel=1e-4)


@pytest.mark.parametrize("S,K,t,r,sig", CASES)
def test_matches_third_diff_of_price(S, K, t, r, sig):
    # Central third difference of the BSM price in sigma.
    h = 1e-2
    c = OptionType.CALL
    third = (price(S, K, t, r, sig + 2 * h, c) - 2 * price(S, K, t, r, sig + h, c)
             + 2 * price(S, K, t, r, sig - h, c)
             - price(S, K, t, r, sig - 2 * h, c)) / (2 * h ** 3)
    assert ultima(S, K, t, r, sig) == pytest.approx(third, rel=3e-2)


def test_call_equals_put():
    assert ultima(100.0, 105.0, 1.0, 0.05, 0.2) == pytest.approx(
        ultima(100.0, 105.0, 1.0, 0.05, 0.2), abs=1e-12)


def test_atm_is_negative():
    # At the money the volga curve is concave in sigma, so ultima < 0.
    assert ultima(100.0, 100.0, 1.0, 0.05, 0.2) < 0.0
