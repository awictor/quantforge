"""Cross-spectral density and magnitude-squared coherence (Welch's method).

Where the power spectrum describes one signal, the *cross*-spectrum describes how two
signals share power at each frequency, and the coherence is the frequency-domain analogue
of a squared correlation: it is 1 where the two are perfectly linearly related at that
frequency and 0 where they are unrelated. Both are estimated by Welch averaging over
overlapping Hann-windowed segments (a single segment gives a meaningless coherence of 1
everywhere, so averaging is essential). Pure standard library on the DFT.
"""

import math

from .spectral import dft


def _welch_spectra(x, y, segment_length, overlap):
    """Averaged auto- and cross-spectra of ``x``, ``y`` over Welch segments.

    Returns ``(freqs, Pxx, Pyy, Pxy)`` where ``Pxy[k]`` is the averaged complex
    cross-spectrum and ``Pxx``/``Pyy`` the averaged auto-spectra, one-sided.
    """
    n = len(x)
    if len(y) != n:
        raise ValueError("x and y must have equal length")
    if n < 16:
        raise ValueError("need at least 16 points")
    L = segment_length if segment_length is not None else max(8, n // 8)
    if L > n:
        raise ValueError("segment_length exceeds the series length")
    if not (0.0 <= overlap < 1.0):
        raise ValueError("overlap must be in [0, 1)")
    step = max(1, int(L * (1.0 - overlap)))
    win = [0.5 - 0.5 * math.cos(2.0 * math.pi * i / (L - 1)) for i in range(L)]
    win_power = sum(w * w for w in win)
    half = L // 2
    pxx = [0.0] * (half + 1)
    pyy = [0.0] * (half + 1)
    pxy = [0j] * (half + 1)
    n_seg = 0
    start = 0
    while start + L <= n:
        xs = x[start:start + L]
        ys = y[start:start + L]
        mx = sum(xs) / L
        my = sum(ys) / L
        fx = dft([(xs[i] - mx) * win[i] for i in range(L)])
        fy = dft([(ys[i] - my) * win[i] for i in range(L)])
        for k in range(half + 1):
            pxx[k] += (abs(fx[k]) ** 2) / win_power
            pyy[k] += (abs(fy[k]) ** 2) / win_power
            pxy[k] += (fx[k] * fy[k].conjugate()) / win_power
        n_seg += 1
        start += step
    if n_seg == 0:
        raise ValueError("no full segments; reduce segment_length")
    freqs = [k / L for k in range(half + 1)]
    pxx = [v / n_seg for v in pxx]
    pyy = [v / n_seg for v in pyy]
    pxy = [v / n_seg for v in pxy]
    return freqs, pxx, pyy, pxy


def cross_spectral_density(x, y, segment_length=None, overlap=0.5):
    """Welch cross-spectral density of two equal-length signals.

    Returns ``(freqs, cross)`` where ``cross[k]`` is the complex averaged cross-spectrum
    at one-sided normalized frequency ``freqs[k]`` in ``[0, 0.5]``. Its magnitude shows
    shared power; its phase, the frequency-dependent lead/lag.
    """
    freqs, _, _, pxy = _welch_spectra(x, y, segment_length, overlap)
    return freqs, pxy


def coherence(x, y, segment_length=None, overlap=0.5):
    """Magnitude-squared coherence ``|Pxy|^2 / (Pxx Pyy)`` in ``[0, 1]`` per frequency.

    The frequency-domain squared correlation: ``1`` where the signals are perfectly
    linearly related at that frequency, ``0`` where unrelated. Returns ``(freqs, coh)``.
    Must average several segments (default ``segment_length = n // 8``) to be meaningful.
    """
    freqs, pxx, pyy, pxy = _welch_spectra(x, y, segment_length, overlap)
    coh = []
    for k in range(len(freqs)):
        denom = pxx[k] * pyy[k]
        coh.append(min(1.0, abs(pxy[k]) ** 2 / denom) if denom > 0 else 0.0)
    return freqs, coh
