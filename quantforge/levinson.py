"""Levinson-Durbin recursion: Toeplitz solves and autoregressive fitting.

A symmetric Toeplitz system -- one whose matrix is constant along each diagonal, as arises
from an autocorrelation sequence -- is solved in ``O(n^2)`` by the Levinson recursion rather
than the ``O(n^3)`` of a general solve. The Durbin specialization fits an autoregressive
model: given the autocorrelation ``r_0..r_p`` it returns the AR coefficients, the prediction
error, and the reflection (PARCOR) coefficients. Pure standard library.
"""


def solve_toeplitz(r, b):
    """Solve ``T x = b`` where ``T`` is the symmetric Toeplitz matrix with first row ``r``.

    ``r[0]`` is the diagonal, ``r[k]`` the ``k``-th off-diagonal (both directions, symmetric).
    ``len(r) == len(b) == n``. ``O(n^2)`` via the Levinson recursion. Raises if a leading
    principal minor is singular (``r[0] == 0`` or a zero prediction error).
    """
    n = len(b)
    if len(r) != n:
        raise ValueError("r and b must have the same length")
    if n == 0:
        return []
    if r[0] == 0:
        raise ValueError("r[0] must be non-zero")
    # Levinson recursion for a symmetric Toeplitz system.
    # y is the Yule-Walker vector solving T_k y = -[r1..rk]; x solves T x = b.
    x = [b[0] / r[0]]
    beta = r[0]
    y = []
    for k in range(1, n):
        # y is the length-k solution of T_k y = -[r1..rk]
        if k == 1:
            alpha = -r[1] / r[0]
            beta = (1 - alpha * alpha) * r[0]
            y = [alpha]
        else:
            num = r[k] + sum(r[k - j] * y[j - 1] for j in range(1, k))
            alpha = -num / beta
            new_y = [y[j] + alpha * y[k - 2 - j] for j in range(k - 1)]
            new_y.append(alpha)
            y = new_y
            beta = (1 - alpha * alpha) * beta
        # update x to length k+1
        mu = (b[k] - sum(r[k - j] * x[j] for j in range(k))) / beta
        new_x = [x[j] + mu * y[k - 1 - j] for j in range(k)]
        new_x.append(mu)
        x = new_x
    return x


def levinson_durbin(autocorr):
    """Fit an AR model from an autocorrelation sequence ``[r_0, r_1, ..., r_p]``.

    Returns ``(ar_coeffs, error, reflection)``: the ``p`` AR coefficients ``a_1..a_p`` such
    that the model predicts ``x_n = sum_k a_k x_{n-k}``, the final prediction error variance,
    and the ``p`` reflection (PARCOR) coefficients. ``O(p^2)``.
    """
    r = list(autocorr)
    p = len(r) - 1
    if p < 0:
        raise ValueError("need at least r_0")
    if r[0] == 0:
        raise ValueError("r_0 must be non-zero")
    a = [0.0] * (p + 1)
    a[0] = 1.0
    err = r[0]
    reflection = []
    for i in range(1, p + 1):
        acc = r[i]
        for j in range(1, i):
            acc += a[j] * r[i - j]
        k = -acc / err
        reflection.append(k)
        new_a = a[:]
        for j in range(1, i):
            new_a[j] = a[j] + k * a[i - j]
        new_a[i] = k
        a = new_a
        err *= (1 - k * k)
        if err <= 0:
            err = 0.0
            break
    # AR prediction coefficients are -a[1..p] (so x_n = sum ar[k] x_{n-1-k})
    ar_coeffs = [-a[j] for j in range(1, p + 1)]
    return ar_coeffs, err, reflection
