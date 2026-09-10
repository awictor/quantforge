"""Craig-Sneyd ADI: second-order-in-time cross-derivative treatment."""

import pytest

from quantforge import adi_two_asset_cs, adi_two_asset, exchange_option


S1, S2, T, R = 100.0, 100.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.25, 0.5


def _pay(s1, s2):
    return max(s1 - s2, 0.0)


@pytest.mark.slow
def test_cs_matches_margrabe():
    mar = exchange_option(S1, S2, T, SIG1, SIG2, RHO)
    cs = adi_two_asset_cs(_pay, S1, S2, T, R, SIG1, SIG2, RHO,
                          n1=140, n2=140, n_time=40)
    assert cs == pytest.approx(mar, abs=1e-2)


@pytest.mark.slow
def test_cs_beats_peaceman_rachford_at_coarse_time():
    # With a mixed derivative, PR is only first-order in time; CS is second
    # order, so at few time steps CS should be markedly more accurate.
    mar = exchange_option(S1, S2, T, SIG1, SIG2, RHO)
    pr = adi_two_asset(_pay, S1, S2, T, R, SIG1, SIG2, RHO,
                       n1=80, n2=80, n_time=10)
    cs = adi_two_asset_cs(_pay, S1, S2, T, R, SIG1, SIG2, RHO,
                          n1=80, n2=80, n_time=10)
    assert abs(cs - mar) < abs(pr - mar)


@pytest.mark.slow
def test_cs_converges_in_space():
    mar = exchange_option(S1, S2, T, SIG1, SIG2, RHO)
    e_lo = abs(adi_two_asset_cs(_pay, S1, S2, T, R, SIG1, SIG2, RHO,
                                n1=60, n2=60, n_time=30) - mar)
    e_hi = abs(adi_two_asset_cs(_pay, S1, S2, T, R, SIG1, SIG2, RHO,
                                n1=140, n2=140, n_time=30) - mar)
    assert e_hi < e_lo


def test_cs_zero_correlation_matches_margrabe():
    mar = exchange_option(S1, S2, T, SIG1, SIG2, 0.0)
    cs = adi_two_asset_cs(_pay, S1, S2, T, R, SIG1, SIG2, 0.0,
                          n1=80, n2=80, n_time=20)
    assert cs == pytest.approx(mar, abs=1e-2)


def test_zero_time_is_payoff():
    assert adi_two_asset_cs(_pay, 110, 100, 0.0, R, SIG1, SIG2, RHO) == \
        pytest.approx(10.0)


def test_bad_prices_raise():
    with pytest.raises(ValueError):
        adi_two_asset_cs(_pay, -1, 100, T, R, SIG1, SIG2, RHO)
