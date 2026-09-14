"""Von Mises distribution: Bessel functions, density, MLE fit."""

import math
import random

import pytest

from quantforge import bessel_i0, bessel_i1, von_mises_pdf, von_mises_fit
from quantforge import resultant_length


def test_bessel_known_values():
    assert abs(bessel_i0(0) - 1.0) < 1e-7
    assert abs(bessel_i0(1) - 1.2660658) < 1e-6
    assert abs(bessel_i0(2) - 2.2795853) < 1e-6
    assert abs(bessel_i1(0)) < 1e-7
    assert abs(bessel_i1(1) - 0.5651591) < 1e-6
    assert abs(bessel_i1(2) - 1.5906369) < 1e-6


def test_pdf_integrates_to_one():
    def integ(mu, kappa, N=20000):
        lo, hi = -math.pi, math.pi
        h = (hi - lo) / N
        s = 0.5 * (von_mises_pdf(lo, mu, kappa) + von_mises_pdf(hi, mu, kappa))
        for i in range(1, N):
            s += von_mises_pdf(lo + i * h, mu, kappa)
        return s * h
    for kappa in [0, 0.5, 2, 5]:
        assert abs(integ(0.3, kappa) - 1.0) < 1e-5


def test_uniform_at_zero_kappa():
    assert abs(von_mises_pdf(1.0, 0, 0) - 1 / (2 * math.pi)) < 1e-9


def test_mode_at_mean():
    mu, kappa = 0.7, 3.0
    peak = von_mises_pdf(mu, mu, kappa)
    for d in [-0.2, -0.1, 0.1, 0.2]:
        assert von_mises_pdf(mu + d, mu, kappa) < peak


def test_fit_recovers_parameters():
    def sample_vm(mu, kappa, n, seed):
        rng = random.Random(seed)
        a = 1 + math.sqrt(1 + 4 * kappa * kappa)
        b = (a - math.sqrt(2 * a)) / (2 * kappa)
        r = (1 + b * b) / (2 * b)
        out = []
        while len(out) < n:
            u1 = rng.random()
            z = math.cos(math.pi * u1)
            f = (1 + r * z) / (r + z)
            c = kappa * (r - f)
            u2 = rng.random()
            if c * (2 - c) - u2 > 0 or math.log(c / u2) + 1 - c >= 0:
                u3 = rng.random()
                theta = mu + (1 if u3 > 0.5 else -1) * math.acos(f)
                out.append(math.atan2(math.sin(theta), math.cos(theta)))
        return out

    data = sample_vm(0.5, 4.0, 5000, 1)
    mu_hat, k_hat = von_mises_fit(data)
    assert abs(mu_hat - 0.5) < 0.1
    assert abs(k_hat - 4.0) < 0.6
    # the fitted kappa satisfies I1/I0 = R
    assert abs(bessel_i1(k_hat) / bessel_i0(k_hat) - resultant_length(data)) < 1e-4


def test_validation():
    with pytest.raises(ValueError):
        von_mises_pdf(0, 0, -1)
    with pytest.raises(ValueError):
        von_mises_fit([1.0])
