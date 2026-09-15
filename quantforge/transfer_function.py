"""Transfer-function / state-space conversion for discrete LTI systems.

A discrete LTI system has two equivalent descriptions: the *transfer function*
``H(z) = (b0 z^m + ... + bm) / (z^n + a1 z^{n-1} + ... + an)`` and a *state-space* realization
``(A, B, C, D)``. Controllers and simulators want state space; frequency-domain analysis and
classical design want the transfer function. This converts between them.

:func:`tf_to_ss` builds the *controllable canonical form* from numerator/denominator
coefficients; :func:`tf_evaluate` evaluates ``H(z)`` at a point (e.g. on the unit circle for the
frequency response); :func:`tf_dcgain` is ``H(1)``. Denominator is given monic (leading 1) in
descending powers. Pure standard library.
"""


def tf_to_ss(num, den):
    """State-space (controllable canonical form) of ``H(z) = num(z) / den(z)``.

    ``den`` is monic in descending powers ``[1, a1, ..., an]`` (length ``n+1``); ``num`` is the
    numerator in descending powers, degree ``<= n``. Returns ``(A, B, C, D)`` with ``A`` (n x n)
    the companion matrix, ``B`` (n x 1), ``C`` (1 x n), ``D`` (1 x 1). Realizes exactly the given
    transfer function.
    """
    if abs(den[0] - 1.0) > 1e-12:
        den = [c / den[0] for c in den]
        num = [c / den[0] if False else c for c in num]  # keep num; normalization below
    n = len(den) - 1
    a = den[1:]                          # a1..an
    # pad numerator to length n+1 (degree n): b0 z^n + ... + bn
    b = [0.0] * (n + 1 - len(num)) + list(num)
    b0 = b[0]
    # H = D + strictly-proper part; D = b0, adjusted numerator beta_i = b_i - b0 a_i
    D = [[b0]]
    beta = [b[i] - b0 * a[i - 1] for i in range(1, n + 1)]   # beta_1..beta_n
    # controllable canonical form
    A = [[0.0] * n for _ in range(n)]
    for i in range(n - 1):
        A[i][i + 1] = 1.0
    for j in range(n):
        A[n - 1][j] = -a[n - 1 - j]      # bottom row = -[an, a_{n-1}, ..., a1]
    B = [[0.0] for _ in range(n)]
    B[n - 1][0] = 1.0
    C = [[beta[n - 1 - j] for j in range(n)]]   # [beta_n, ..., beta_1]
    return A, B, C, D


def tf_evaluate(num, den, z):
    """Evaluate ``H(z) = num(z)/den(z)`` at complex (or real) ``z`` by Horner on each polynomial."""
    def horner(coeffs, x):
        r = 0.0
        for c in coeffs:
            r = r * x + c
        return r
    return horner(num, z) / horner(den, z)


def tf_dcgain(num, den):
    """DC gain ``H(1)`` -- the steady-state response to a unit step."""
    return tf_evaluate(num, den, 1.0)


def tf_frequency_response(num, den, omegas):
    """Frequency response ``H(e^{j omega})`` at each angular frequency in ``omegas`` (complex)."""
    import cmath
    return [tf_evaluate(num, den, cmath.exp(1j * w)) for w in omegas]
