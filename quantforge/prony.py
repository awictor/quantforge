"""Prony's method: fit a sum of complex exponentials to uniformly-sampled data.

Prony's method (1795) recovers the parameters of a signal modelled as a sum of ``p`` damped
complex exponentials,

    y_k = sum_{i=1}^{p} a_i * z_i**k ,   k = 0, 1, ..., N-1 ,

from equally-spaced samples. It is exact when the data really is such a sum (up to noise) and is
the classical route to fitting decaying sinusoids, extracting modal frequencies/damping from a
ring-down, and rational-function partial-fraction recovery. The trick: the exponentials satisfy
a linear recurrence whose characteristic polynomial has the ``z_i`` as roots, so you solve a
linear system for the recurrence coefficients, root the polynomial to get the ``z_i``, then a
linear least-squares for the amplitudes ``a_i``.

Returns the modes ``z_i`` (complex) and amplitudes ``a_i`` (complex). For a real sampling
interval ``dt`` the continuous rate is ``log(z_i) / dt`` (real part = damping, imag = angular
frequency). Pure standard library on top of :func:`quantforge.polyroots.polynomial_roots`.
"""

from .polyroots import polynomial_roots


def _solve_complex(A, b):
    # Gaussian elimination with partial pivoting over complex numbers
    n = len(A)
    M = [[complex(A[i][j]) for j in range(n)] + [complex(b[i])] for i in range(n)]
    for k in range(n):
        piv = max(range(k, n), key=lambda i: abs(M[i][k]))
        if abs(M[piv][k]) < 1e-300:
            raise ValueError("singular system in Prony solve")
        M[k], M[piv] = M[piv], M[k]
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]
    x = [0j] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = s / M[i][i]
    return x


def _lstsq_complex(A, b):
    # least squares via the normal equations A^H A x = A^H b
    m = len(A)
    n = len(A[0])
    AhA = [[sum(A[k][i].conjugate() * A[k][j] for k in range(m)) for j in range(n)]
           for i in range(n)]
    Ahb = [sum(A[k][i].conjugate() * b[k] for k in range(m)) for i in range(n)]
    return _solve_complex(AhA, Ahb)


def prony(y, p):
    """Fit ``p`` complex exponentials to the sample sequence ``y`` by Prony's method.

    Requires ``len(y) >= 2 p``. Returns ``(modes, amplitudes)`` where ``y_k ~ sum amplitudes[i] *
    modes[i]**k``. Modes and amplitudes are complex; a real signal yields conjugate pairs.
    """
    N = len(y)
    if p < 1:
        raise ValueError("need p >= 1")
    if N < 2 * p:
        raise ValueError("need at least 2p samples")
    yc = [complex(v) for v in y]

    # Step 1: linear-prediction coefficients. Solve the Hankel system
    #   sum_{j=0}^{p-1} c_j y_{k+j} = -y_{k+p}   for k = 0 .. N-p-1  (least squares).
    A = [[yc[k + j] for j in range(p)] for k in range(N - p)]
    b = [-yc[k + p] for k in range(N - p)]
    c = _lstsq_complex(A, b)

    # Step 2: characteristic polynomial z^p + c_{p-1} z^{p-1} + ... + c_0, root for modes.
    # polynomial_roots wants highest-degree-first coefficients.
    coeffs = [1.0] + [c[p - 1 - i] for i in range(p)]
    modes = polynomial_roots(coeffs)

    # Step 3: amplitudes by least squares on the Vandermonde system y_k = sum a_i modes_i^k.
    V = [[modes[i] ** k for i in range(p)] for k in range(N)]
    amps = _lstsq_complex(V, yc)
    return modes, amps


def prony_reconstruct(modes, amps, n):
    """Reconstruct ``n`` samples ``y_k = sum amps[i] * modes[i]**k`` from Prony parameters."""
    return [sum(amps[i] * modes[i] ** k for i in range(len(modes))) for k in range(n)]
