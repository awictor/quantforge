"""Fast linear recurrences via matrix exponentiation.

A linear recurrence ``x_n = c_1 x_{n-1} + ... + c_k x_{n-k}`` advances one step by a fixed
``k x k`` companion matrix, so the ``n``-th term is that matrix raised to the ``n``-th power
applied to the initial state. Exponentiation by squaring computes the power in
``O(k^3 log n)`` -- far faster than iterating ``n`` steps once ``n`` is large -- and works
over the integers, the reals, or modulo an integer. Pure standard library.
"""


def matrix_power(matrix, power, mod=None):
    """Raise a square ``matrix`` to a non-negative integer ``power`` by squaring.

    ``power == 0`` returns the identity. With ``mod`` set, every entry is reduced modulo it
    (integer matrices only). ``O(k^3 log power)``.
    """
    if power < 0:
        raise ValueError("power must be non-negative")
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError("matrix must be square")

    def mul(a, b):
        out = [[0] * n for _ in range(n)]
        for i in range(n):
            ai = a[i]
            oi = out[i]
            for t in range(n):
                v = ai[t]
                if v == 0:
                    continue
                bt = b[t]
                for j in range(n):
                    oi[j] += v * bt[j]
            if mod is not None:
                for j in range(n):
                    oi[j] %= mod
        return out

    # identity
    result = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    base = [row[:] for row in matrix]
    if mod is not None:
        base = [[x % mod for x in row] for row in base]
        result = [[x % mod for x in row] for row in result]
    p = power
    while p > 0:
        if p & 1:
            result = mul(result, base)
        p >>= 1
        if p:
            base = mul(base, base)
    return result


def linear_recurrence_nth(coeffs, initial, n, mod=None):
    """The ``n``-th term (0-indexed) of ``x_i = sum_j coeffs[j] * x_{i-1-j}``.

    ``coeffs`` are the ``k`` recurrence coefficients (``coeffs[0]`` multiplies the most
    recent term); ``initial`` are the first ``k`` terms ``x_0 .. x_{k-1}``. For ``n < k`` the
    initial term is returned directly. With ``mod`` set, the result is reduced modulo it.
    ``O(k^3 log n)``.
    """
    k = len(coeffs)
    if len(initial) != k:
        raise ValueError("initial must have the same length as coeffs")
    if n < 0:
        raise ValueError("n must be non-negative")
    if k == 0:
        raise ValueError("need at least one coefficient")
    if n < k:
        return initial[n] % mod if mod is not None else initial[n]
    # companion matrix: top row = coeffs, sub-diagonal identity
    comp = [[0] * k for _ in range(k)]
    for j in range(k):
        comp[0][j] = coeffs[j]
    for i in range(1, k):
        comp[i][i - 1] = 1
    # state vector [x_{k-1}, x_{k-2}, ..., x_0]; comp^(n-k+1) advances it to x_n on top
    powered = matrix_power(comp, n - (k - 1), mod=mod)
    state = list(reversed(initial))    # [x_{k-1}, ..., x_0]
    val = sum(powered[0][j] * state[j] for j in range(k))
    return val % mod if mod is not None else val


def fibonacci(n, mod=None):
    """The ``n``-th Fibonacci number (``F_0 = 0, F_1 = 1``) via matrix power.

    Optional modular reduction. ``O(log n)`` matrix multiplies.
    """
    if n < 0:
        raise ValueError("n must be non-negative")
    return linear_recurrence_nth([1, 1], [0, 1], n, mod=mod)
