"""MUSIC: MUltiple SIgnal Classification frequency estimation.

Prony (:mod:`quantforge.prony`) and the matrix pencil (:mod:`quantforge.matrix_pencil`) recover
signal modes directly. MUSIC (Schmidt 1986) instead builds a *pseudospectrum* whose peaks mark
the frequencies of ``p`` sinusoids in noise. It splits the data's subspace into a ``p``-dimensional
signal subspace and its orthogonal complement (the noise subspace) via an SVD of the Hankel data
matrix; a sinusoid's steering vector is orthogonal to the noise subspace, so the reciprocal of
its projection onto that subspace spikes at the true frequencies. MUSIC resolves closely-spaced
tones far below the Fourier (1/N) limit, which is why it is a staple of array processing and
radar.

:func:`music_pseudospectrum` evaluates the pseudospectrum on a frequency grid;
:func:`music_frequencies` returns the ``p`` peak frequencies. Real-valued input. Pure standard
library on top of :mod:`quantforge.svd`.
"""

import cmath
import math

from .svd import svd


def _hankel(y, L):
    rows = len(y) - L + 1
    return [[y[i + j] for j in range(L)] for i in range(rows)]


def _noise_subspace(y, p, L):
    # SVD of the real Hankel matrix; right-singular vectors beyond the top p span the
    # noise subspace (columns of V for the smallest singular values).
    H = _hankel([float(v) for v in y], L)
    U, S, V = svd(H)          # V is L x L, columns = right-singular vectors
    n = len(S)
    # order indices by descending singular value; keep the tail (noise) columns
    order = sorted(range(n), key=lambda t: -S[t])
    noise_cols = order[p:]
    # noise subspace basis vectors (each length L)
    return [[V[i][t] for i in range(L)] for t in noise_cols]


def music_pseudospectrum(y, p, freqs, L=None):
    """MUSIC pseudospectrum of ``y`` at each angular frequency in ``freqs`` (radians/sample).

    ``p`` is the number of (real) sinusoids -- use ``2 * (number of tones)`` since each real tone
    is a conjugate pair. ``L`` is the Hankel window length (default ``len(y)//2``). Returns a list
    of pseudospectrum values ``1 / ||E_n^H a(omega)||^2``; peaks sit at the tone frequencies.
    """
    n = len(y)
    if L is None:
        L = n // 2
    if not (p < L <= n):
        raise ValueError("require p < L <= len(y)")
    En = _noise_subspace(y, p, L)      # list of noise basis vectors, each length L
    out = []
    for w in freqs:
        # steering vector a(w) = [e^{-i w k}], k = 0..L-1
        a = [cmath.exp(-1j * w * k) for k in range(L)]
        # projection energy onto noise subspace = sum_v |v^H a|^2  (v real)
        energy = 0.0
        for v in En:
            proj = sum(v[k] * a[k] for k in range(L))
            energy += abs(proj) ** 2
        out.append(1.0 / energy if energy > 1e-300 else float("inf"))
    return out


def music_frequencies(y, p, L=None, grid=4000):
    """Estimate the ``p`` dominant frequencies of ``y`` as the pseudospectrum peaks.

    Scans ``[0, pi]`` on ``grid`` points and returns the frequencies of the ``ceil(p/2)`` largest
    local maxima (each real tone shows one peak in ``[0, pi]``). Sorted ascending.
    """
    n = len(y)
    if L is None:
        L = n // 2
    freqs = [math.pi * i / grid for i in range(grid + 1)]
    ps = music_pseudospectrum(y, p, freqs, L)
    # local maxima
    peaks = []
    for i in range(1, len(ps) - 1):
        if ps[i] > ps[i - 1] and ps[i] >= ps[i + 1]:
            peaks.append((ps[i], freqs[i]))
    peaks.sort(key=lambda t: -t[0])
    n_tones = (p + 1) // 2
    chosen = sorted(f for _, f in peaks[:n_tones])
    return chosen
