"""Savitzky-Golay smoothing and differentiation filter.

Fits a low-degree polynomial by least squares to a sliding window of points and
evaluates it (or its derivative) at the window center. This smooths noise while
preserving peak shape far better than a moving average, and exactly reproduces any
polynomial up to the fit degree. The convolution coefficients are obtained from the
normal equations of the window's Vandermonde matrix; the derivative uses the
polynomial's derivative coefficient. Pure standard library on top of the matrix
inverse.
"""

from .portopt import _invert


def savgol_coeffs(window, degree, deriv=0):
    """Savitzky-Golay convolution coefficients for a window of odd length ``window``.

    Returns the ``window`` weights that, dotted with the windowed samples, give the
    fitted value (``deriv = 0``) or the ``deriv``-th derivative at the window center.
    ``window`` must be odd and larger than ``degree``. For unit spacing; scale a
    derivative by ``1 / h^deriv`` for spacing ``h``.
    """
    if window % 2 == 0 or window < 1:
        raise ValueError("window must be a positive odd integer")
    if degree >= window:
        raise ValueError("degree must be less than the window length")
    if deriv > degree:
        raise ValueError("deriv must not exceed degree")
    half = window // 2
    # Vandermonde A[i][j] = z_i^j, z_i = -half..half.
    A = [[float(z) ** j for j in range(degree + 1)]
         for z in range(-half, half + 1)]
    # Normal equations: (A' A)^{-1} A'. The center-fit coefficients are row `deriv`
    # of that pseudo-inverse, times deriv! (the derivative of z^deriv at z=0).
    At = [[A[i][j] for i in range(window)] for j in range(degree + 1)]
    AtA = [[sum(At[i][k] * A[k][j] for k in range(window)) for j in range(degree + 1)]
           for i in range(degree + 1)]
    inv = _invert(AtA)
    # (A'A)^{-1} A' has shape (degree+1) x window; row `deriv` gives the coefficients.
    import math
    fact = math.factorial(deriv)
    coeffs = []
    for i in range(window):
        c = sum(inv[deriv][k] * At[k][i] for k in range(degree + 1))
        coeffs.append(c * fact)
    return coeffs


def savgol_filter(data, window, degree, deriv=0):
    """Apply a Savitzky-Golay filter to ``data``.

    Smooths (``deriv = 0``) or differentiates the series with a length-``window``,
    degree-``degree`` polynomial fit. Interior points use the centered convolution;
    the ``half`` points at each end are fitted from the nearest full window (so the
    output has the same length as the input). Exactly reproduces polynomials up to
    ``degree``.
    """
    n = len(data)
    half = window // 2
    if n < window:
        raise ValueError("data shorter than the window")
    coeffs = savgol_coeffs(window, degree, deriv)
    out = [0.0] * n
    for t in range(half, n - half):
        out[t] = sum(coeffs[k] * data[t - half + k] for k in range(window))
    # Edges: fit the polynomial on the first/last full window and evaluate off-center.
    for end, base in ((0, 0), (1, n - window)):
        for offset in range(half):
            pos = offset if end == 0 else n - half + offset
            # Coefficients for evaluating at z = (pos - center_of_window).
            center = base + half
            z = pos - center
            c = _offset_coeffs(window, degree, deriv, z)
            out[pos] = sum(c[k] * data[base + k] for k in range(window))
    return out


def _offset_coeffs(window, degree, deriv, z):
    """Savitzky-Golay coefficients evaluated at offset ``z`` from the window center."""
    import math
    half = window // 2
    A = [[float(zz) ** j for j in range(degree + 1)]
         for zz in range(-half, half + 1)]
    At = [[A[i][j] for i in range(window)] for j in range(degree + 1)]
    AtA = [[sum(At[i][k] * A[k][j] for k in range(window)) for j in range(degree + 1)]
           for i in range(degree + 1)]
    inv = _invert(AtA)
    # Evaluate the deriv-th derivative of the fitted polynomial at z.
    coeffs = [0.0] * window
    for p in range(deriv, degree + 1):
        # d^deriv/dz^deriv z^p = p!/(p-deriv)! z^(p-deriv)
        fac = math.factorial(p) / math.factorial(p - deriv)
        basis = fac * (z ** (p - deriv))
        for i in range(window):
            coeffs[i] += basis * sum(inv[p][k] * At[k][i] for k in range(degree + 1))
    return coeffs
