"""Tests for Breeden-Litzenberger risk-neutral density extraction.

We build a Black-Scholes call curve (whose risk-neutral law is lognormal),
recover the density/CDF from it, and check that the density integrates to one,
reprices vanilla payoffs, and that the CDF has the right tail limits.
"""

import math

import pytest

from quantforge import (
    call_price, put_price,
    risk_neutral_density, risk_neutral_cdf, price_from_density, density_total_mass,
)


def _bsm_call_curve(S=100, t=1.0, r=0.05, sigma=0.2, lo=1.0, hi=300.0, step=0.5):
    strikes = [lo + i * step for i in range(int((hi - lo) / step) + 1)]
    calls = [call_price(S, K, t, r, sigma) for K in strikes]
    return strikes, calls


def test_density_integrates_to_one():
    strikes, calls = _bsm_call_curve()
    assert density_total_mass(strikes, calls, 1.0, 0.05) == pytest.approx(1.0, abs=1e-3)


def test_density_is_nonnegative():
    strikes, calls = _bsm_call_curve()
    _, dens = risk_neutral_density(strikes, calls, 1.0, 0.05)
    assert all(d >= 0 for d in dens)


def test_density_matches_lognormal_pdf():
    # The BL density of a BSM curve is the lognormal risk-neutral pdf.
    S, t, r, sigma = 100, 1.0, 0.05, 0.2
    strikes, calls = _bsm_call_curve(S, t, r, sigma, step=0.25)
    mids, dens = risk_neutral_density(strikes, calls, t, r)

    def lognormal_pdf(x):
        # Risk-neutral density of S_t: lognormal with mean (r-0.5 sig^2)t.
        m = math.log(S) + (r - 0.5 * sigma * sigma) * t
        s = sigma * math.sqrt(t)
        return math.exp(-((math.log(x) - m) ** 2) / (2 * s * s)) / (x * s * math.sqrt(2 * math.pi))

    for k, d in zip(mids, dens):
        if 60 <= k <= 160:  # where the density has meaningful mass
            assert d == pytest.approx(lognormal_pdf(k), abs=2e-4)


def test_reprices_vanilla_put():
    S, t, r, sigma = 100, 1.0, 0.05, 0.2
    strikes, calls = _bsm_call_curve(S, t, r, sigma)
    val = price_from_density(strikes, calls, t, r, lambda ST: max(110 - ST, 0.0))
    assert val == pytest.approx(put_price(S, 110, t, r, sigma), abs=1e-2)


def test_reprices_vanilla_call():
    S, t, r, sigma = 100, 1.0, 0.05, 0.2
    strikes, calls = _bsm_call_curve(S, t, r, sigma)
    val = price_from_density(strikes, calls, t, r, lambda ST: max(ST - 95, 0.0))
    assert val == pytest.approx(call_price(S, 95, t, r, sigma), abs=1e-2)


def test_reprices_digital():
    # A cash-or-nothing call paying 1 if S_T > 100.
    S, t, r, sigma = 100, 1.0, 0.05, 0.2
    strikes, calls = _bsm_call_curve(S, t, r, sigma)
    val = price_from_density(strikes, calls, t, r, lambda ST: 1.0 if ST > 100 else 0.0)
    # Closed form: e^{-rt} N(d2).
    from quantforge.mathfns import norm_cdf
    d2 = (math.log(S / 100) + (r - 0.5 * sigma * sigma) * t) / (sigma * math.sqrt(t))
    expected = math.exp(-r * t) * norm_cdf(d2)
    assert val == pytest.approx(expected, abs=5e-3)


def test_cdf_tail_limits():
    strikes, calls = _bsm_call_curve()
    mids, cdf = risk_neutral_cdf(strikes, calls, 1.0, 0.05)
    # Deep OTM strike -> CDF near 0; deep ITM -> near 1.
    assert cdf[0] < 0.01
    assert cdf[-1] > 0.99
    # Monotone non-decreasing.
    assert all(cdf[i] <= cdf[i + 1] + 1e-9 for i in range(len(cdf) - 1))


def test_requires_three_strikes():
    with pytest.raises(ValueError):
        risk_neutral_density([90, 100], [12, 6], 1.0, 0.05)
