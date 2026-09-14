"""Sample-rate conversion: sinc interpolation and anti-aliased up/downsampling.

Changing a signal's sample rate cannot be done by naively dropping or repeating
samples -- downsampling that way folds high frequencies back as aliases, and upsampling
that way leaves spectral images. The correct operations low-pass filter at the new
Nyquist limit. ``upsample`` inserts zeros and interpolates with an anti-imaging filter;
``downsample`` anti-alias filters then keeps every ``factor``-th sample;
``resample_rational`` combines them for a rational ratio. ``sinc_interp`` reconstructs
the underlying band-limited continuous signal at arbitrary points (Whittaker-Shannon).
Pure standard library on top of the windowed-sinc FIR design.
"""

import math

from .fir_filter import fir_lowpass, fir_apply


def sinc_interp(x, positions):
    """Whittaker-Shannon reconstruction of a band-limited signal at arbitrary ``positions``.

    Treats ``x`` as samples at integer indices ``0..len(x)-1`` of a signal band-limited
    to the Nyquist frequency and evaluates ``sum_n x[n] sinc(t - n)`` at each ``t`` in
    ``positions``. Exact at integer positions; the ideal interpolation between them.
    """
    n = len(x)
    if n == 0:
        raise ValueError("input must be non-empty")
    out = []
    for t in positions:
        acc = 0.0
        for k in range(n):
            d = t - k
            if d == 0.0:
                acc += x[k]
            else:
                pi_d = math.pi * d
                acc += x[k] * math.sin(pi_d) / pi_d
        out.append(acc)
    return out


def upsample(x, factor, numtaps=65):
    """Upsample ``x`` by an integer ``factor`` with an anti-imaging low-pass filter.

    Inserts ``factor - 1`` zeros between samples and low-pass filters at the original
    Nyquist (cutoff ``1/factor``), scaling by ``factor`` to preserve amplitude. Returns a
    signal ``factor`` times as long. ``numtaps`` sizes the interpolation filter (odd).
    """
    if factor < 1:
        raise ValueError("factor must be >= 1")
    n = len(x)
    if n == 0:
        raise ValueError("input must be non-empty")
    if factor == 1:
        return list(x)
    expanded = [0.0] * (n * factor)
    for i in range(n):
        expanded[i * factor] = x[i]
    taps = fir_lowpass(numtaps, 1.0 / factor)
    filtered = fir_apply(taps, expanded)
    return [v * factor for v in filtered]


def downsample(x, factor, numtaps=65):
    """Downsample ``x`` by an integer ``factor`` with an anti-alias low-pass filter.

    Low-pass filters at the new Nyquist (cutoff ``1/factor``) to prevent aliasing, then
    keeps every ``factor``-th sample. Returns a signal ``ceil(len(x)/factor)`` long.
    ``numtaps`` sizes the anti-alias filter (odd).
    """
    if factor < 1:
        raise ValueError("factor must be >= 1")
    n = len(x)
    if n == 0:
        raise ValueError("input must be non-empty")
    if factor == 1:
        return list(x)
    taps = fir_lowpass(numtaps, 1.0 / factor)
    filtered = fir_apply(taps, x)
    return [filtered[i] for i in range(0, n, factor)]


def resample_rational(x, up, down, numtaps=65):
    """Resample ``x`` by the rational ratio ``up/down`` (upsample then downsample).

    Interpolates by ``up`` and decimates by ``down`` through a single anti-alias/anti-
    imaging low-pass at cutoff ``1/max(up, down)``, so both imaging and aliasing are
    suppressed. Returns a signal about ``len(x) * up / down`` long.
    """
    if up < 1 or down < 1:
        raise ValueError("up and down must be >= 1")
    n = len(x)
    if n == 0:
        raise ValueError("input must be non-empty")
    if up == 1 and down == 1:
        return list(x)
    expanded = [0.0] * (n * up)
    for i in range(n):
        expanded[i * up] = x[i]
    taps = fir_lowpass(numtaps, 1.0 / max(up, down))
    filtered = fir_apply(taps, expanded)
    scaled = [v * up for v in filtered]
    return [scaled[i] for i in range(0, len(scaled), down)]
