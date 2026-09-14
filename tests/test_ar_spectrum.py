"""Parametric (AR) power-spectral-density estimation: Burg + Yule-Walker."""

import math
import random

import pytest

from quantforge import ar_psd, burg, ar_spectrum


def _ar2_series(phi, N, seed):
    rng = random.Random(seed)
    x = [0.0, 0.0]
    for t in range(2, N):
        x.append(phi[0] * x[t - 1] + phi[1] * x[t - 2] + rng.gauss(0.0, 1.0))
    return x[100:]      # drop burn-in


def test_burg_recovers_ar2_coefficients():
    phi = [0.75, -0.5]
    x = _ar2_series(phi, 2000, seed=7)
    m = burg(x, 2)
    assert abs(m["coefficients"][0] - phi[0]) < 0.05
    assert abs(m["coefficients"][1] - phi[1]) < 0.05


def test_reflection_coefficients_stable():
    x = _ar2_series([0.75, -0.5], 2000, seed=11)
    m = burg(x, 4)
    assert all(abs(k) < 1.0 for k in m["reflection"])   # stable model


def test_spectrum_peaks_at_sinusoids():
    rng = random.Random(3)
    f1, f2, M = 0.1, 0.25, 128
    sig = [math.sin(2 * math.pi * f1 * n) + math.sin(2 * math.pi * f2 * n) +
           0.1 * rng.gauss(0.0, 1.0) for n in range(M)]
    freqs = [i / 500 for i in range(251)]
    psd = ar_spectrum(sig, 20, freqs, method="burg")
    # both tones show up as local maxima (the second is weaker but still a peak)
    peaks = [freqs[i] for i in range(1, len(psd) - 1)
             if psd[i] > psd[i - 1] and psd[i] > psd[i + 1] and psd[i] > max(psd) * 1e-3]
    assert any(abs(p - f1) < 0.01 for p in peaks)
    assert any(abs(p - f2) < 0.01 for p in peaks)


def test_ar1_psd_shape():
    freqs = [i / 500 for i in range(251)]
    # positive phi -> low-pass (more power at DC than Nyquist)
    lp = ar_psd([0.9], 1.0, freqs)
    assert lp[0] > lp[-1]
    # negative phi -> high-pass
    hp = ar_psd([-0.9], 1.0, freqs)
    assert hp[-1] > hp[0]


def test_yule_walker_method_matches_burg_on_long_series():
    x = _ar2_series([0.6, -0.3], 3000, seed=5)
    freqs = [i / 100 for i in range(51)]
    pb = ar_spectrum(x, 2, freqs, method="burg")
    py = ar_spectrum(x, 2, freqs, method="yule_walker")
    # both estimate the same spectrum on a long series
    peak_b = max(range(len(pb)), key=lambda i: pb[i])
    peak_y = max(range(len(py)), key=lambda i: py[i])
    assert abs(peak_b - peak_y) <= 2


def test_validation():
    with pytest.raises(ValueError):
        burg([1.0, 2.0, 3.0], 0)               # order < 1
    with pytest.raises(ValueError):
        burg([1.0, 2.0], 5)                    # order >= len
    with pytest.raises(ValueError):
        ar_spectrum([1.0, 2.0, 3.0, 4.0], 2, [0.1], method="nope")
    with pytest.raises(ValueError):
        ar_psd([0.5], -1.0, [0.1])             # negative noise variance
