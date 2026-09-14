"""Hilbert transform and the analytic signal (FFT method).

The analytic signal of a real series is a complex series whose real part is the
original and whose imaginary part is its Hilbert transform -- a 90-degree phase shift
of every frequency component. It turns an oscillation into a rotating phasor, from
which the *envelope* (instantaneous amplitude) and the *instantaneous phase and
frequency* fall out directly. This is the standard demodulation tool for AM/FM
signals, seismic and vibration analysis, and empirical-mode decomposition.

The transform is computed the usual way: take the FFT, zero the negative frequencies
and double the positive ones (Marple's one-sided construction), then invert. Requires
a power-of-two length. Pure standard library (wraps :mod:`quantforge.fft`).
"""

import cmath
import math

from .fft import fft, ifft


def analytic_signal(x):
    """Analytic signal of a real sequence ``x`` (length a power of two).

    Returns a list of complex values ``z`` with ``z.real == x`` (to floating error)
    and ``z.imag`` the Hilbert transform of ``x``. The one-sided spectrum is formed by
    keeping the DC and Nyquist bins, doubling the positive-frequency bins, and zeroing
    the negative-frequency bins, then inverting the FFT.
    """
    n = len(x)
    if n == 0 or (n & (n - 1)) != 0:
        raise ValueError("length must be a positive power of two")
    X = fft(x)
    # Multiplier h: 1 at DC and Nyquist, 2 on positive frequencies, 0 on negatives.
    h = [0.0] * n
    if n % 2 == 0:
        h[0] = 1.0
        h[n // 2] = 1.0
        for k in range(1, n // 2):
            h[k] = 2.0
    else:  # unreachable for power-of-two n>=2, kept for completeness
        h[0] = 1.0
        for k in range(1, (n + 1) // 2):
            h[k] = 2.0
    return ifft([X[k] * h[k] for k in range(n)])


def hilbert_transform(x):
    """Discrete Hilbert transform of ``x``: the imaginary part of its analytic signal.

    A 90-degree phase shift of every frequency component (cosine -> sine). Length must
    be a power of two.
    """
    return [z.imag for z in analytic_signal(x)]


def envelope(x):
    """Instantaneous-amplitude envelope ``|analytic_signal(x)|``.

    For an amplitude-modulated carrier this recovers the modulating envelope; for a
    pure sinusoid it is a constant equal to the amplitude. Length a power of two.
    """
    return [abs(z) for z in analytic_signal(x)]


def instantaneous_phase(x):
    """Unwrapped instantaneous phase (radians) of the analytic signal of ``x``.

    ``atan2(imag, real)`` per sample, unwrapped so it is continuous rather than
    jumping by ``2*pi``. Length a power of two.
    """
    z = analytic_signal(x)
    raw = [math.atan2(v.imag, v.real) for v in z]
    out = [raw[0]] if raw else []
    offset = 0.0
    for i in range(1, len(raw)):
        d = raw[i] - raw[i - 1]
        while d > math.pi:
            d -= 2.0 * math.pi
        while d < -math.pi:
            d += 2.0 * math.pi
        offset = out[i - 1] + d
        out.append(offset)
    return out


def instantaneous_frequency(x, dt=1.0):
    """Instantaneous frequency (cycles per unit time) from the analytic signal.

    The derivative of the unwrapped phase, ``d(phase)/dt / (2*pi)``, by central
    differences on the interior and one-sided differences at the ends. For a pure tone
    this is flat at the tone's frequency. ``dt`` is the sample spacing; length a power
    of two.
    """
    phase = instantaneous_phase(x)
    n = len(phase)
    if n < 2:
        raise ValueError("need at least two samples")
    freq = [0.0] * n
    two_pi = 2.0 * math.pi
    freq[0] = (phase[1] - phase[0]) / dt / two_pi
    freq[-1] = (phase[-1] - phase[-2]) / dt / two_pi
    for i in range(1, n - 1):
        freq[i] = (phase[i + 1] - phase[i - 1]) / (2.0 * dt) / two_pi
    return freq
