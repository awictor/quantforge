"""Heston pathwise delta vs a Fourier finite-difference."""

import pytest

from quantforge import (
    OptionType,
    heston_pathwise_delta,
    heston_price,
)


S, K, T, R = 100.0, 100.0, 1.0, 0.03
V0, KAPPA, THETA, XI, RHO = 0.04, 1.5, 0.04, 0.5, -0.7


def _fourier_delta(ot):
    h = 0.5
    up = heston_price(S + h, K, T, R, V0, KAPPA, THETA, XI, RHO, ot)
    dn = heston_price(S - h, K, T, R, V0, KAPPA, THETA, XI, RHO, ot)
    return (up - dn) / (2 * h)


@pytest.mark.slow
def test_call_delta_matches_fourier_fd():
    pw = heston_pathwise_delta(S, K, T, R, V0, KAPPA, THETA, XI, RHO,
                               OptionType.CALL, n_steps=100, n_paths=120_000,
                               seed=1)
    assert pw.price == pytest.approx(_fourier_delta(OptionType.CALL),
                                     abs=3.0 * pw.std_error + 5e-3)


@pytest.mark.slow
def test_put_delta_matches_fourier_fd():
    pw = heston_pathwise_delta(S, K, T, R, V0, KAPPA, THETA, XI, RHO,
                               OptionType.PUT, n_steps=100, n_paths=120_000,
                               seed=2)
    assert pw.price == pytest.approx(_fourier_delta(OptionType.PUT),
                                     abs=3.0 * pw.std_error + 5e-3)


def test_call_delta_in_unit_interval():
    pw = heston_pathwise_delta(S, K, T, R, V0, KAPPA, THETA, XI, RHO,
                               OptionType.CALL, n_steps=50, n_paths=20_000,
                               seed=3)
    assert 0.0 < pw.price < 1.0


def test_put_delta_negative():
    pw = heston_pathwise_delta(S, K, T, R, V0, KAPPA, THETA, XI, RHO,
                               OptionType.PUT, n_steps=50, n_paths=20_000,
                               seed=4)
    assert pw.price < 0.0


def test_bad_params_raise():
    with pytest.raises(ValueError):
        heston_pathwise_delta(-1, K, T, R, V0, KAPPA, THETA, XI, RHO)
