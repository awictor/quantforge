"""Special functions: incomplete gamma and beta, digamma, inverse erf."""

import math

import pytest

from quantforge import gammainc, gammaincc, betainc, digamma, erfinv

EULER_GAMMA = 0.5772156649015329


def test_gamma_p_plus_q_is_one():
    for a, x in [(0.5, 0.3), (2.0, 1.5), (5.0, 12.0), (3.0, 0.1), (10.0, 8.0)]:
        assert abs(gammainc(a, x) + gammaincc(a, x) - 1.0) < 1e-14


def test_gamma_half_is_erf():
    # P(1/2, x) = erf(sqrt(x))
    for x in (0.25, 1.0, 4.0, 9.0):
        assert abs(gammainc(0.5, x) - math.erf(math.sqrt(x))) < 1e-12


def test_gamma_one_is_exponential_cdf():
    # P(1, x) = 1 - e^{-x}
    for x in (0.5, 2.0, 5.0):
        assert abs(gammainc(1.0, x) - (1.0 - math.exp(-x))) < 1e-12


def test_gamma_limits():
    assert gammainc(3.0, 0.0) == 0.0
    assert gammaincc(3.0, 0.0) == 1.0
    assert gammainc(2.0, 100.0) > 1.0 - 1e-12


def test_beta_known_values():
    assert abs(betainc(2.0, 2.0, 0.5) - 0.5) < 1e-12          # symmetric a=b
    assert abs(betainc(1.0, 1.0, 0.3) - 0.3) < 1e-12          # I_x(1,1) = x
    assert betainc(2.0, 3.0, 0.0) == 0.0
    assert betainc(2.0, 3.0, 1.0) == 1.0


def test_beta_symmetry():
    # I_x(a, b) = 1 - I_{1-x}(b, a)
    for a, b, x in [(2.0, 3.0, 0.4), (0.5, 2.5, 0.7), (4.0, 1.5, 0.2)]:
        assert abs(betainc(a, b, x) - (1.0 - betainc(b, a, 1.0 - x))) < 1e-12


def test_beta_monotone_in_x():
    prev = -1.0
    for k in range(11):
        x = k / 10.0
        v = betainc(2.0, 5.0, x)
        assert v >= prev
        prev = v


def test_digamma_special_values():
    assert abs(digamma(1.0) + EULER_GAMMA) < 1e-9
    assert abs(digamma(0.5) - (-EULER_GAMMA - 2.0 * math.log(2.0))) < 1e-9
    # Recurrence psi(x+1) = psi(x) + 1/x
    for x in (0.7, 2.3, 5.0):
        assert abs(digamma(x + 1.0) - digamma(x) - 1.0 / x) < 1e-10


def test_erfinv_roundtrip():
    for y in (-0.9, -0.3, 0.1, 0.5, 0.9, 0.99):
        assert abs(math.erf(erfinv(y)) - y) < 1e-14
    assert erfinv(0.0) == 0.0


def test_erfinv_matches_probit_scale():
    # erfinv(y) = ppf((1+y)/2) / sqrt(2); check against a known normal quantile.
    from quantforge.mathfns import norm_ppf
    y = 0.6
    assert abs(erfinv(y) - norm_ppf((1 + y) / 2) / math.sqrt(2)) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        gammainc(0.0, 1.0)
    with pytest.raises(ValueError):
        gammainc(1.0, -1.0)
    with pytest.raises(ValueError):
        betainc(1.0, 1.0, 1.5)
    with pytest.raises(ValueError):
        betainc(-1.0, 1.0, 0.5)
    with pytest.raises(ValueError):
        digamma(0.0)
    with pytest.raises(ValueError):
        erfinv(1.5)
