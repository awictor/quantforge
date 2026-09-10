"""Carr-Madan FFT strip: one transform prices the whole smile."""

import bisect
import cmath
import math

import pytest

from quantforge import (
    OptionType,
    carr_madan_strip,
    carr_madan_smile_strip,
    cgmy_price,
    nig_price,
)
from quantforge.carrmadan import _fft
from quantforge.cgmy import _cgmy_psi
from quantforge.nig import _nig_psi


def test_fft_matches_naive_dft():
    x = [1.0, 2.0, -1.0, 3.0, 0.5, -2.0, 1.5, 0.0]
    n = len(x)
    dft = [sum(x[k] * cmath.exp(-2j * math.pi * j * k / n) for k in range(n))
           for j in range(n)]
    got = _fft(x)
    assert max(abs(a - b) for a, b in zip(got, dft)) < 1e-12


def test_fft_inverse_roundtrip():
    x = [1.0, 2.0, -1.0, 3.0, 0.5, -2.0, 1.5, 0.0]
    back = _fft(_fft(x), inverse=True)
    assert max(abs(a - b.real) for a, b in zip(x, back)) < 1e-12


def test_fft_length_must_be_power_of_two():
    with pytest.raises(ValueError):
        _fft([1.0, 2.0, 3.0])


def test_cgmy_strip_matches_per_strike():
    S, t, r, q = 100.0, 1.0, 0.03, 0.0
    C, G, M, Y = 4.0, 5.0, 10.0, 0.5
    strikes, calls = carr_madan_strip(S, t, r, q,
                                      lambda u: _cgmy_psi(u, C, G, M, Y))
    # Grid strikes near the money should match the per-strike GL pricer closely.
    i = bisect.bisect_left(strikes, 100.0)
    for j in (i - 1, i, i + 1):
        K = strikes[j]
        exact = cgmy_price(S, K, t, r, C, G, M, Y, OptionType.CALL, q=q)
        assert calls[j] == pytest.approx(exact, abs=1e-4)


def test_nig_strip_matches_per_strike():
    S, t, r, q = 100.0, 1.0, 0.03, 0.0
    a, b, d = 15.0, -5.0, 0.5
    strikes, calls = carr_madan_strip(S, t, r, q,
                                      lambda u: _nig_psi(u, a, b, d))
    i = bisect.bisect_left(strikes, 100.0)
    for j in (i - 1, i, i + 1):
        K = strikes[j]
        exact = nig_price(S, K, t, r, a, b, d, OptionType.CALL, q=q)
        assert calls[j] == pytest.approx(exact, abs=1e-4)


def test_smile_strip_returns_sorted_window():
    S, t, r, q = 100.0, 0.5, 0.03, 0.0
    sm = carr_madan_smile_strip(S, t, r, q,
                                lambda u: _nig_psi(u, 15.0, -6.0, 0.4),
                                k_lo=-0.2, k_hi=0.2)
    lms = [lm for lm, _ in sm]
    assert lms == sorted(lms)
    assert all(-0.2 <= lm <= 0.2 for lm in lms)
    assert all(iv > 0 for _, iv in sm)
    # beta < 0 => downward skew across the window.
    assert sm[0][1] > sm[-1][1]


def test_strip_requires_power_of_two():
    with pytest.raises(ValueError):
        carr_madan_strip(100.0, 1.0, 0.03, 0.0,
                         lambda u: _nig_psi(u, 15.0, -5.0, 0.5), n_fft=1000)
