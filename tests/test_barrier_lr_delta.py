"""Barrier-option delta by the likelihood-ratio method (barrier_lr_delta)."""

import pytest

from quantforge import barrier_lr_delta, barrier_mc, OptionType


S, K, T, R, SIG, NS = 100.0, 100.0, 1.0, 0.05, 0.2, 50


def _crn_fd(H, ot, barrier, h=0.5, seed=100):
    # Common-random-number finite difference of the discretely-monitored price.
    up = barrier_mc(S + h, K, H, T, R, SIG, ot, barrier, n_steps=NS,
                    n_paths=200_000, seed=seed, brownian_bridge=False)
    dn = barrier_mc(S - h, K, H, T, R, SIG, ot, barrier, n_steps=NS,
                    n_paths=200_000, seed=seed, brownian_bridge=False)
    return (up.price - dn.price) / (2 * h)


@pytest.mark.slow
def test_down_out_call_matches_crn_fd():
    lr = barrier_lr_delta(S, K, 90.0, T, R, SIG, OptionType.CALL, "down-out",
                          n_steps=NS, n_paths=200_000, seed=1)
    fd = _crn_fd(90.0, OptionType.CALL, "down-out", seed=1)
    assert lr.price == pytest.approx(fd, abs=3.0 * lr.std_error + 5e-3)


@pytest.mark.slow
def test_up_in_put_matches_crn_fd():
    lr = barrier_lr_delta(S, K, 110.0, T, R, SIG, OptionType.PUT, "up-in",
                          n_steps=NS, n_paths=200_000, seed=2)
    fd = _crn_fd(110.0, OptionType.PUT, "up-in", seed=2)
    assert lr.price == pytest.approx(fd, abs=3.0 * lr.std_error + 5e-3)


def test_down_out_call_delta_positive():
    lr = barrier_lr_delta(S, K, 90.0, T, R, SIG, OptionType.CALL, "down-out",
                          n_steps=NS, n_paths=60_000, seed=3)
    assert lr.price > 0.0


def test_bad_barrier_raises():
    with pytest.raises(ValueError):
        barrier_lr_delta(S, K, 90.0, T, R, SIG, barrier="sideways-out")


def test_bad_H_raises():
    with pytest.raises(ValueError):
        barrier_lr_delta(S, K, -1.0, T, R, SIG)
