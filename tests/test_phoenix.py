"""Phoenix autocallable note (Monte Carlo)."""

import pytest

from quantforge import phoenix_autocall_mc as ph


S, T, R, SIG = 100, 3.0, 0.03, 0.25
OBS = [1.0, 2.0, 3.0]
AC, CB, CPN, PB = 110, 80, 0.08, 70


def test_price_positive():
    assert ph(S, T, R, SIG, OBS, AC, CB, CPN, protection_barrier=PB,
              n_paths=40000, seed=1) > 0


def test_memory_at_least_no_memory():
    pm = ph(S, T, R, SIG, OBS, AC, CB, CPN, protection_barrier=PB, memory=True,
            n_paths=40000, seed=1)
    pn = ph(S, T, R, SIG, OBS, AC, CB, CPN, protection_barrier=PB, memory=False,
            n_paths=40000, seed=1)
    assert pm >= pn


def test_lower_coupon_barrier_pays_more():
    lo = ph(S, T, R, SIG, OBS, AC, 60, CPN, protection_barrier=PB, n_paths=40000, seed=2)
    hi = ph(S, T, R, SIG, OBS, AC, 95, CPN, protection_barrier=PB, n_paths=40000, seed=2)
    assert lo > hi


def test_higher_coupon_higher_value():
    assert (ph(S, T, R, SIG, OBS, AC, CB, 0.15, protection_barrier=PB, n_paths=40000, seed=3)
            > ph(S, T, R, SIG, OBS, AC, CB, 0.05, protection_barrier=PB, n_paths=40000, seed=3))


def test_no_protection_at_least_with_protection():
    with_pb = ph(S, T, R, SIG, OBS, AC, CB, CPN, protection_barrier=PB, n_paths=40000, seed=1)
    no_pb = ph(S, T, R, SIG, OBS, AC, CB, CPN, protection_barrier=None, n_paths=40000, seed=1)
    assert no_pb >= with_pb


def test_validation():
    with pytest.raises(ValueError):
        ph(-1, T, R, SIG, OBS, AC, CB, CPN)
    with pytest.raises(ValueError):
        ph(S, T, R, SIG, [], AC, CB, CPN)
