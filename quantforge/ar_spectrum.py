"""Parametric (autoregressive) power-spectral-density estimation.

The periodogram estimates a spectrum by transforming the data directly; a *parametric*
estimate instead fits an autoregressive model and reads its spectrum off the
coefficients. For a short record with a few sharp spectral peaks -- sinusoids in noise,
a resonance, a formant -- the AR spectrum resolves them far better than a periodogram of
the same length, because it does not assume the signal is zero outside the window. Two
routes to the coefficients are provided: Yule-Walker (via the existing fit) and Burg's
method, which minimizes forward and backward prediction error and is the standard choice
for short series. Pure standard library.
"""

import math

from .ar_model import fit_ar_yule_walker


def ar_psd(coeffs, noise_variance, freqs):
    """AR power spectral density ``sigma^2 / |1 - sum phi_k e^{-i 2 pi f k}|^2``.

    ``freqs`` are normalized frequencies in cycles/sample (``0`` to ``0.5`` is DC to
    Nyquist). ``coeffs`` are the AR coefficients ``phi_1..phi_p`` and
    ``noise_variance`` the innovation variance. Returns the PSD at each frequency.
    """
    if noise_variance < 0.0:
        raise ValueError("noise_variance must be non-negative")
    p = len(coeffs)
    out = []
    for f in freqs:
        # Denominator A(f) = 1 - sum phi_k e^{-i 2 pi f k}.
        re = 1.0
        im = 0.0
        for k in range(1, p + 1):
            ang = -2.0 * math.pi * f * k
            re -= coeffs[k - 1] * math.cos(ang)
            im -= coeffs[k - 1] * math.sin(ang)
        denom = re * re + im * im
        out.append(noise_variance / denom if denom > 0.0 else float("inf"))
    return out


def burg(x, order):
    """Fit AR coefficients by Burg's method (minimizes forward+backward error).

    Returns ``{"coefficients": [phi_1..phi_p], "noise_variance": sigma2,
    "reflection": [k_1..k_p]}``. Burg estimates the reflection coefficients directly
    from the data (never forming autocovariances), which gives sharper, more stable
    spectra than Yule-Walker on short records. ``order`` must be ``>= 1`` and less than
    ``len(x)``.
    """
    n = len(x)
    if order < 1:
        raise ValueError("order must be >= 1")
    if n <= order:
        raise ValueError("series too short for the requested order")
    f = [float(v) for v in x]          # forward prediction errors
    b = [float(v) for v in x]          # backward prediction errors
    a = [1.0]                          # AR polynomial coefficients (a[0] == 1)
    # Total power (den) drives the reflection-coefficient denominator.
    dk = sum(2.0 * v * v for v in x) - f[0] * f[0] - b[n - 1] * b[n - 1]
    var = sum(v * v for v in x) / n
    reflection = []
    for m in range(order):
        # Reflection coefficient k = -2 * sum f[i] b[i-1] / dk.
        num = 0.0
        for i in range(m + 1, n):
            num += f[i] * b[i - 1]
        k = -2.0 * num / dk if dk != 0.0 else 0.0
        reflection.append(-k)          # report the usual-sign reflection coefficient
        # Update the AR polynomial: a_new[i] = a[i] + k * a[m+1-i], with a[m+1] == 0.
        a_prev = a + [0.0]              # index m+1 now valid and zero
        a = a_prev[:]
        for i in range(1, m + 2):
            a[i] = a_prev[i] + k * a_prev[m + 1 - i]
        var *= (1.0 - k * k)
        # Update forward/backward errors in place (walk high to low to avoid clobber).
        new_f = f[:]
        new_b = b[:]
        for i in range(n - 1, m, -1):
            new_f[i] = f[i] + k * b[i - 1]
            new_b[i] = b[i - 1] + k * f[i]
        f, b = new_f, new_b
        dk = (1.0 - k * k) * dk - f[m + 1] * f[m + 1] - b[n - 1] * b[n - 1]
    coeffs = [-a[i] for i in range(1, order + 1)]   # x_t = sum phi_i x_{t-i} + eps
    return {"coefficients": coeffs, "noise_variance": var, "reflection": reflection}


def ar_spectrum(x, order, freqs, method="burg"):
    """Parametric PSD of ``x`` from an AR(``order``) fit, evaluated at ``freqs``.

    ``method`` is ``"burg"`` (default, best for short records) or ``"yule_walker"``.
    Returns the PSD at each normalized frequency in ``freqs``.
    """
    if method == "burg":
        m = burg(x, order)
    elif method == "yule_walker":
        m = fit_ar_yule_walker(x, order)
    else:
        raise ValueError("method must be 'burg' or 'yule_walker'")
    return ar_psd(m["coefficients"], m["noise_variance"], freqs)
