"""Bivariate copulas: Gaussian, Clayton, Gumbel."""

import pytest

from quantforge import (
    gaussian_copula, clayton_copula, gumbel_copula,
    clayton_lower_tail_dependence, gumbel_upper_tail_dependence,
    clayton_theta_from_tau, gumbel_theta_from_tau,
    frank_copula, gaussian_copula_joint_default, first_to_default_probability,
    vasicek_loss_cdf, vasicek_loss_quantile, cdo_tranche_expected_loss,
)


def test_boundary_zero():
    assert gaussian_copula(0.5, 0, 0.5) == 0
    assert clayton_copula(0.5, 0, 2) == 0
    assert gumbel_copula(0.5, 0, 2) == 0


def test_boundary_uniform_margin():
    assert gaussian_copula(0.5, 1, 0.5) == pytest.approx(0.5)
    assert clayton_copula(0.5, 1, 2) == pytest.approx(0.5)
    assert gumbel_copula(0.5, 1, 2) == pytest.approx(0.5)


def test_independence_limits():
    assert gaussian_copula(0.5, 0.6, 0.0) == pytest.approx(0.3, abs=1e-9)
    assert clayton_copula(0.5, 0.6, 1e-8) == pytest.approx(0.3, abs=1e-6)
    assert gumbel_copula(0.5, 0.6, 1.0) == pytest.approx(0.3, abs=1e-12)


def test_frechet_bounds_and_positive_dependence():
    u, v = 0.4, 0.7
    lo, hi = max(u + v - 1, 0), min(u, v)
    for val in (gaussian_copula(u, v, 0.5), clayton_copula(u, v, 2),
                gumbel_copula(u, v, 2)):
        assert lo - 1e-9 <= val <= hi + 1e-9
        assert val > u * v   # positive dependence


def test_clayton_lower_tail_dependence():
    assert clayton_lower_tail_dependence(2) == pytest.approx(2 ** (-0.5))
    assert 0 < clayton_lower_tail_dependence(2) < 1


def test_gumbel_upper_tail_dependence():
    assert gumbel_upper_tail_dependence(2) == pytest.approx(2 - 2 ** 0.5)
    assert gumbel_upper_tail_dependence(1.0) == pytest.approx(0.0)


def test_theta_from_tau_round_trips():
    th = clayton_theta_from_tau(0.5)
    assert th / (th + 2) == pytest.approx(0.5)
    thg = gumbel_theta_from_tau(0.5)
    assert 1 - 1 / thg == pytest.approx(0.5)


def test_vasicek_cdf_monotone_and_bounded():
    xs = [0.01, 0.05, 0.1, 0.2, 0.5, 0.9]
    cdfs = [vasicek_loss_cdf(x, 0.05, 0.2) for x in xs]
    assert all(cdfs[i] < cdfs[i + 1] for i in range(len(cdfs) - 1))
    assert all(0 <= c <= 1 for c in cdfs)


def test_vasicek_quantile_inverts_and_monotone():
    L = vasicek_loss_quantile(0.99, 0.05, 0.2)
    assert vasicek_loss_cdf(L, 0.05, 0.2) == pytest.approx(0.99, abs=1e-6)
    assert vasicek_loss_quantile(0.999, 0.05, 0.2) > vasicek_loss_quantile(0.99, 0.05, 0.2)
    assert vasicek_loss_quantile(0.99, 0.10, 0.2) > vasicek_loss_quantile(0.99, 0.05, 0.2)
    assert vasicek_loss_quantile(0.99, 0.05, 0.4) > vasicek_loss_quantile(0.99, 0.05, 0.2)


def test_cdo_full_capital_structure_mean_is_pd():
    assert cdo_tranche_expected_loss(0, 1, 0.05, 0.2) == pytest.approx(0.05, abs=1e-3)


def test_cdo_tranche_seniority_ordering():
    eq = cdo_tranche_expected_loss(0, 0.03, 0.05, 0.2)
    mez = cdo_tranche_expected_loss(0.03, 0.07, 0.05, 0.2)
    sen = cdo_tranche_expected_loss(0.07, 0.15, 0.05, 0.2)
    assert eq > mez > sen
    assert all(0 <= x <= 1 for x in (eq, mez, sen))


def test_vasicek_cdo_validation():
    with pytest.raises(ValueError):
        vasicek_loss_cdf(0.5, 1.5, 0.2)
    with pytest.raises(ValueError):
        cdo_tranche_expected_loss(0.5, 0.3, 0.05, 0.2)


def test_frank_boundary_and_independence():
    assert frank_copula(0.5, 0, 2) == 0
    assert frank_copula(0.5, 1, 2) == pytest.approx(0.5)
    assert frank_copula(0.4, 0.7, 1e-9) == pytest.approx(0.28, abs=1e-6)


def test_frank_symmetric_and_signed_dependence():
    assert frank_copula(0.4, 0.7, 3) == pytest.approx(frank_copula(0.7, 0.4, 3))
    assert frank_copula(0.4, 0.7, 5) > 0.28   # positive dependence
    assert frank_copula(0.4, 0.7, -5) < 0.28  # negative dependence


def test_joint_default_above_product_when_correlated():
    pd1, pd2 = 0.05, 0.08
    assert gaussian_copula_joint_default(pd1, pd2, 0.5) >= pd1 * pd2
    assert gaussian_copula_joint_default(pd1, pd2, 0.0) == pytest.approx(pd1 * pd2, abs=1e-9)


def test_first_to_default_bounds_and_monotonicity():
    pd1, pd2 = 0.05, 0.08
    f = first_to_default_probability(pd1, pd2, 0.5)
    assert max(pd1, pd2) <= f <= min(pd1 + pd2, 1)
    assert first_to_default_probability(pd1, pd2, 0.8) < \
        first_to_default_probability(pd1, pd2, 0.0)
    assert first_to_default_probability(pd1, pd2, 0.0) == pytest.approx(
        pd1 + pd2 - pd1 * pd2, abs=1e-9)


def test_default_copula_validation():
    with pytest.raises(ValueError):
        gaussian_copula_joint_default(1.5, 0.5, 0.3)


def test_validation():
    with pytest.raises(ValueError):
        clayton_copula(0.5, 0.5, -1)
    with pytest.raises(ValueError):
        gumbel_copula(0.5, 0.5, 0.5)
    with pytest.raises(ValueError):
        clayton_theta_from_tau(1.0)
