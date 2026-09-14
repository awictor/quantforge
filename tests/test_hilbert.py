"""Hilbert transform, analytic signal, envelope, instantaneous phase/frequency."""

import math

import pytest

from quantforge import (
    analytic_signal,
    hilbert_transform,
    envelope,
    instantaneous_phase,
    instantaneous_frequency,
)


def _cos(f, N):
    return [math.cos(2 * math.pi * f * n / N) for n in range(N)]


def test_analytic_real_part_is_input():
    N = 256
    x = _cos(3, N)
    z = analytic_signal(x)
    assert max(abs(z[n].real - x[n]) for n in range(N)) < 1e-12


def test_hilbert_of_cosine_is_sine():
    N = 256
    x = _cos(3, N)
    z = analytic_signal(x)
    sin_ref = [math.sin(2 * math.pi * 3 * n / N) for n in range(N)]
    # interior avoids circular-convolution edge wrap
    assert max(abs(z[n].imag - sin_ref[n]) for n in range(20, N - 20)) < 1e-10


def test_hilbert_transform_matches_imag_part():
    N = 128
    x = _cos(5, N)
    ht = hilbert_transform(x)
    z = analytic_signal(x)
    assert max(abs(ht[n] - z[n].imag) for n in range(N)) < 1e-12


def test_envelope_of_pure_tone_is_amplitude():
    N = 256
    x = [2.0 * math.cos(2 * math.pi * 5 * n / N) for n in range(N)]
    env = envelope(x)
    mid = env[50:200]
    assert abs(sum(mid) / len(mid) - 2.0) < 1e-6


def test_envelope_recovers_am_modulation():
    N = 256
    fc, fm = 40, 2
    am = [(1.0 + 0.5 * math.cos(2 * math.pi * fm * n / N)) *
          math.cos(2 * math.pi * fc * n / N) for n in range(N)]
    env = envelope(am)
    mod = [1.0 + 0.5 * math.cos(2 * math.pi * fm * n / N) for n in range(N)]
    assert max(abs(env[n] - mod[n]) for n in range(30, N - 30)) < 1e-3


def test_instantaneous_frequency_of_tone_is_flat():
    N = 256
    f = 3
    x = _cos(f, N)
    inst = instantaneous_frequency(x, dt=1.0)
    mid = inst[40:200]
    assert abs(sum(mid) / len(mid) - f / N) < 1e-6
    assert max(mid) - min(mid) < 1e-6          # flat


def test_phase_is_unwrapped_monotone_for_tone():
    N = 256
    x = _cos(3, N)
    ph = instantaneous_phase(x)
    # a positive-frequency tone has monotonically increasing unwrapped phase (interior)
    assert all(ph[n + 1] > ph[n] for n in range(30, N - 30))


def test_validation():
    with pytest.raises(ValueError):
        analytic_signal([1.0, 2.0, 3.0])       # not a power of two
    with pytest.raises(ValueError):
        analytic_signal([])
