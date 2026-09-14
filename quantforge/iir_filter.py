"""Butterworth IIR filter design (bilinear transform, biquad cascade).

An IIR filter reaches a given frequency selectivity with far fewer coefficients than a
FIR filter by feeding output back into the recursion -- the trade is nonlinear phase.
The Butterworth response is maximally flat in the passband (no ripple), which makes it
the default general-purpose choice. These design it the textbook way: place the analog
prototype poles on a circle, warp the cutoff for the bilinear transform, map each analog
pole pair to a digital second-order section (biquad), and cascade the sections. Working
in second-order sections rather than one high-order polynomial keeps the recursion
numerically stable. Pure standard library.
"""

import math


def _prototype_poles(order):
    """Left-half-plane poles of the order-``order`` Butterworth analog prototype (unit cutoff).

    Returns a list of complex poles; real (unpaired) pole included when ``order`` is odd.
    """
    poles = []
    for k in range(order):
        theta = math.pi * (2 * k + 1) / (2 * order) + math.pi / 2.0
        poles.append(complex(math.cos(theta), math.sin(theta)))
    return poles


def _biquad_from_analog_lowpass(p_pair, wc, highpass=False):
    """Digital biquad (b, a) from an analog low-pass pole pair via the bilinear transform.

    ``p_pair`` is one or two prototype poles (a real pole is handled as a first-order
    section padded to second order). ``wc`` is the pre-warped cutoff. Set ``highpass`` to
    map ``s -> wc / s`` first (low-pass prototype to high-pass).
    """
    # Scale prototype (unit cutoff) to wc; for high-pass use the s -> wc/s transform,
    # which sends each pole p to wc/p and adds zeros at the origin.
    if highpass:
        scaled = [wc / p for p in p_pair]
    else:
        scaled = [wc * p for p in p_pair]
    # Bilinear transform s = 2*fs*(z-1)/(z+1) with fs = 1 (pre-warping already applied),
    # i.e. s = K*(1-z^-1)/(1+z^-1), K = 2. Build the section numerator/denominator.
    K = 2.0
    # Denominator from poles: prod (K*(1-z^-1) - scaled_pole*(1+z^-1)).
    den = [1.0 + 0j]
    for p in scaled:
        # factor in z^-1: (K - p) + (-K - p) z^-1
        factor = [(K - p), (-K - p)]
        den = _poly_mul(den, factor)
    # Numerator: low-pass has zeros at s = inf -> z = -1, i.e. (1 + z^-1) per pole;
    # high-pass has zeros at s = 0 -> z = +1, i.e. (1 - z^-1) per pole.
    num = [1.0 + 0j]
    zero_factor = [1.0, -1.0] if highpass else [1.0, 1.0]
    for _ in scaled:
        num = _poly_mul(num, zero_factor)
    # Normalize so a[0] == 1.
    a0 = den[0]
    b = [(c / a0).real for c in num]
    a = [(c / a0).real for c in den]
    return b, a


def _poly_mul(p, q):
    out = [0j] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            out[i + j] += pi * qj
    return out


def _design(order, cutoff, highpass):
    if order < 1:
        raise ValueError("order must be >= 1")
    if not (0.0 < cutoff < 0.5):
        raise ValueError("cutoff must be in (0, 0.5) cycles/sample")
    # Pre-warp the cutoff for the bilinear transform (fs = 1).
    wc = 2.0 * math.tan(math.pi * cutoff)
    poles = _prototype_poles(order)
    # Pair conjugate poles into biquads; a lone real pole (odd order) goes solo.
    sections = []
    used = [False] * order
    for i in range(order):
        if used[i]:
            continue
        pi = poles[i]
        if abs(pi.imag) < 1e-12:
            sections.append(_biquad_from_analog_lowpass([pi], wc, highpass))
            used[i] = True
        else:
            # find conjugate
            for j in range(i + 1, order):
                if not used[j] and abs(poles[j] - pi.conjugate()) < 1e-9:
                    sections.append(_biquad_from_analog_lowpass([pi, poles[j]], wc, highpass))
                    used[i] = used[j] = True
                    break
    # Normalize overall gain to unity at the reference frequency (DC for LP, Nyquist for HP).
    ref = 1.0 if not highpass else -1.0     # z = 1 (DC) or z = -1 (Nyquist)
    gain = 1.0
    for b, a in sections:
        num = sum(b[k] * ref ** k for k in range(len(b)))
        den = sum(a[k] * ref ** k for k in range(len(a)))
        gain *= num / den
    if sections and gain != 0.0:
        b0, a0 = sections[0]
        sections[0] = ([c / gain for c in b0], a0)
    return sections


def butter_lowpass(order, cutoff):
    """Design a Butterworth low-pass as a cascade of biquad sections.

    ``cutoff`` is the -3 dB frequency in cycles/sample, in ``(0, 0.5)``. Returns a list
    of ``(b, a)`` second-order sections (each a length-3 numerator and denominator with
    ``a[0] == 1``), with the overall DC gain normalized to 1.
    """
    return _design(order, cutoff, highpass=False)


def butter_highpass(order, cutoff):
    """Design a Butterworth high-pass as a cascade of biquad sections.

    ``cutoff`` is the -3 dB frequency in cycles/sample, in ``(0, 0.5)``. Returns a list
    of ``(b, a)`` second-order sections with the Nyquist gain normalized to 1.
    """
    return _design(order, cutoff, highpass=True)


def sosfilt(sections, x):
    """Filter ``x`` through a cascade of biquad ``(b, a)`` sections (Direct Form II transposed).

    Runs the signal through each second-order section in turn. Returns a list the same
    length as ``x``.
    """
    y = [float(v) for v in x]
    for b, a in sections:
        y = _biquad_apply(b, a, y)
    return y


def _biquad_apply(b, a, x):
    # Direct Form II transposed, second order (pad to length 3).
    b0, b1, b2 = (b + [0.0, 0.0, 0.0])[:3]
    a0, a1, a2 = (a + [0.0, 0.0, 0.0])[:3]
    out = []
    z1 = 0.0
    z2 = 0.0
    for xn in x:
        yn = b0 * xn + z1
        z1 = b1 * xn - a1 * yn + z2
        z2 = b2 * xn - a2 * yn
        out.append(yn)
    return out


def iir_frequency_response(sections, freqs):
    """Magnitude response ``|H(f)|`` of a biquad cascade at normalized ``freqs``.

    ``freqs`` in cycles/sample (``0`` to ``0.5``). Returns the product of the section
    magnitudes at each frequency -- useful for verifying a design's passband/stopband.
    """
    out = []
    for f in freqs:
        z_inv = complex(math.cos(-2 * math.pi * f), math.sin(-2 * math.pi * f))
        h = 1.0 + 0j
        for b, a in sections:
            num = sum((b + [0.0, 0.0, 0.0])[k] * z_inv ** k for k in range(3))
            den = sum((a + [0.0, 0.0, 0.0])[k] * z_inv ** k for k in range(3))
            h *= num / den
        out.append(abs(h))
    return out
