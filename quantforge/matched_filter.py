"""Matched filtering and peak detection: find a known template in noisy data.

The matched filter is the linear filter that maximizes signal-to-noise ratio when
detecting a *known* waveform buried in white noise -- it correlates the data against a
time-reversed copy of the template. Its output peaks where the template best aligns, so
it is the optimal detector for radar/sonar returns, gravitational-wave templates, a known
spike shape in a sensor trace, or any "find this pattern" task. ``matched_filter``
produces the detection statistic; ``normalized_matched_filter`` scales it to a
correlation coefficient; ``find_peaks`` and ``detect_template`` locate the hits. Pure
standard library.
"""

import math


def matched_filter(x, template):
    """Matched-filter response: correlate ``x`` against ``template`` at every offset.

    Returns a list of length ``len(x) - len(template) + 1`` whose entry ``i`` is the dot
    product of ``template`` with the window ``x[i:i+len(template)]``. The maximum locates
    the offset where the template best aligns -- the maximum-SNR detector for a known
    shape in white noise.
    """
    n = len(x)
    m = len(template)
    if m == 0 or n == 0:
        raise ValueError("inputs must be non-empty")
    if m > n:
        raise ValueError("template longer than signal")
    out = []
    for i in range(n - m + 1):
        out.append(sum(x[i + j] * template[j] for j in range(m)))
    return out


def normalized_matched_filter(x, template):
    """Normalized matched filter: the correlation coefficient at each offset, in ``[-1, 1]``.

    Removes each window's mean and the template mean, then divides by the norms, so the
    response is a dimensionless correlation independent of signal amplitude and offset --
    ``1`` where the window is a scaled, shifted copy of the template. Returns a list of
    length ``len(x) - len(template) + 1``.
    """
    n = len(x)
    m = len(template)
    if m == 0 or n == 0:
        raise ValueError("inputs must be non-empty")
    if m > n:
        raise ValueError("template longer than signal")
    tmean = sum(template) / m
    tc = [t - tmean for t in template]
    tnorm = math.sqrt(sum(v * v for v in tc))
    if tnorm == 0.0:
        raise ValueError("template has zero variance")
    out = []
    for i in range(n - m + 1):
        win = x[i:i + m]
        wmean = sum(win) / m
        wc = [w - wmean for w in win]
        wnorm = math.sqrt(sum(v * v for v in wc))
        if wnorm == 0.0:
            out.append(0.0)
        else:
            out.append(sum(wc[j] * tc[j] for j in range(m)) / (wnorm * tnorm))
    return out


def find_peaks(x, height=None, distance=1):
    """Indices of local maxima in ``x``, optionally filtered by height and spacing.

    A point is a peak if it is strictly greater than both neighbors. ``height`` drops
    peaks below that threshold; ``distance`` enforces a minimum gap between reported peaks
    by keeping the taller peak when two fall within ``distance`` samples. Returns a sorted
    list of indices.
    """
    n = len(x)
    if distance < 1:
        raise ValueError("distance must be >= 1")
    candidates = [i for i in range(1, n - 1) if x[i] > x[i - 1] and x[i] > x[i + 1]]
    if height is not None:
        candidates = [i for i in candidates if x[i] >= height]
    if distance == 1 or not candidates:
        return candidates
    # Greedily keep the tallest peaks, suppressing any within `distance`.
    order = sorted(candidates, key=lambda i: x[i], reverse=True)
    kept = []
    for i in order:
        if all(abs(i - j) >= distance for j in kept):
            kept.append(i)
    return sorted(kept)


def detect_template(x, template, threshold=0.8, distance=None):
    """Offsets where ``template`` occurs in ``x`` (normalized correlation above ``threshold``).

    Runs :func:`normalized_matched_filter` and returns the peak offsets whose correlation
    exceeds ``threshold`` (a coefficient in ``[-1, 1]``), separated by at least
    ``distance`` samples (default ``len(template)``, so overlapping detections of the same
    hit collapse to one). Each returned offset is the start index of a match.

    Because the correlation is amplitude-normalized, short or featureless templates match
    many noise windows by shape alone; use a longer, distinctive template (and a higher
    ``threshold``) when false positives matter.
    """
    if distance is None:
        distance = len(template)
    corr = normalized_matched_filter(x, template)
    return find_peaks(corr, height=threshold, distance=max(1, distance))
