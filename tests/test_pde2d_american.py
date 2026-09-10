"""Two-asset American ADI (early-exercise projection) vs closed forms."""

import pytest

from quantforge import (
    adi_two_asset_american,
    best_of_call_closed,
    worst_of_call_closed,
)


S1, S2, K, T, R = 100.0, 100.0, 100.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.25, 0.3


def _best(s1, s2):
    return max(max(s1, s2) - K, 0.0)


def _worst(s1, s2):
    return max(min(s1, s2) - K, 0.0)


@pytest.mark.slow
def test_european_best_of_matches_closed_form():
    eu = adi_two_asset_american(_best, S1, S2, T, R, SIG1, SIG2, RHO,
                                n1=100, n2=100, n_time=60, american=False)
    cf = best_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    assert eu == pytest.approx(cf, abs=8e-2)


def test_european_worst_of_converges_to_closed_form():
    # Exact Stulz closed form (deterministic) as the PDE reference -- was
    # previously the Monte Carlo worst_of_call, whose noise occasionally exceeded
    # the 3e-2 band.
    cf = worst_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    eu = adi_two_asset_american(_worst, S1, S2, T, R, SIG1, SIG2, RHO,
                                n1=100, n2=100, n_time=100, american=False)
    assert eu == pytest.approx(cf, abs=3e-2)


def test_american_at_least_european():
    eu = adi_two_asset_american(_best, S1, S2, T, R, SIG1, SIG2, RHO,
                                n1=60, n2=60, n_time=40, american=False)
    am = adi_two_asset_american(_best, S1, S2, T, R, SIG1, SIG2, RHO,
                                n1=60, n2=60, n_time=40, american=True)
    assert am >= eu - 1e-6


def test_dividends_give_early_exercise_premium():
    eu = adi_two_asset_american(_best, S1, S2, T, R, SIG1, SIG2, RHO,
                                q1=0.06, q2=0.06, n1=70, n2=70, n_time=50,
                                american=False)
    am = adi_two_asset_american(_best, S1, S2, T, R, SIG1, SIG2, RHO,
                                q1=0.06, q2=0.06, n1=70, n2=70, n_time=50,
                                american=True)
    assert am > eu


def test_american_floored_at_intrinsic():
    # Deep in the money, an American best-of is worth at least its intrinsic.
    v = adi_two_asset_american(_best, 150, 140, T, R, SIG1, SIG2, RHO,
                               q1=0.06, q2=0.06, n1=60, n2=60, n_time=40,
                               american=True)
    assert v >= max(150, 140) - K - 1e-6


def test_zero_time_is_payoff():
    assert adi_two_asset_american(_best, 130, 120, 0.0, R, SIG1, SIG2, RHO) == \
        pytest.approx(30.0)
