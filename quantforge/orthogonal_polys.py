"""Classical orthogonal polynomials via their three-term recurrences.

The Hermite, Laguerre, and Chebyshev families are the orthogonal polynomials attached to the
Gaussian, exponential, and arcsine weights. They are the eigenfunctions of the quantum
harmonic oscillator (Hermite), the radial hydrogen problem (Laguerre), and the workhorses of
minimax approximation and spectral methods (Chebyshev). Each is evaluated here by its stable
three-term recurrence -- no factorials, no cancellation -- so degrees into the hundreds stay
accurate. Pure standard library.
"""

import math


def hermite_h(n, x):
    """Physicists' Hermite polynomial ``H_n(x)`` (weight ``e^{-x^2}`` on the real line).

    ``H_0 = 1``, ``H_1 = 2x``, ``H_{k+1} = 2x H_k - 2k H_{k-1}``. ``H_2 = 4x^2 - 2``.
    """
    if n < 0:
        raise ValueError("hermite_h requires n >= 0")
    if n == 0:
        return 1.0
    h_prev = 1.0
    h_cur = 2.0 * x
    for k in range(1, n):
        h_next = 2.0 * x * h_cur - 2.0 * k * h_prev
        h_prev, h_cur = h_cur, h_next
    return h_cur


def hermite_he(n, x):
    """Probabilists' Hermite polynomial ``He_n(x)`` (weight ``e^{-x^2/2}``).

    ``He_0 = 1``, ``He_1 = x``, ``He_{k+1} = x He_k - k He_{k-1}``. Related to the physicists'
    form by ``He_n(x) = 2^{-n/2} H_n(x / sqrt 2)``.
    """
    if n < 0:
        raise ValueError("hermite_he requires n >= 0")
    if n == 0:
        return 1.0
    h_prev = 1.0
    h_cur = x
    for k in range(1, n):
        h_next = x * h_cur - k * h_prev
        h_prev, h_cur = h_cur, h_next
    return h_cur


def laguerre_l(n, x, alpha=0.0):
    """Generalized Laguerre polynomial ``L_n^{(alpha)}(x)`` (weight ``x^alpha e^{-x}``).

    ``L_0 = 1``, ``L_1 = 1 + alpha - x``, and
    ``(k+1) L_{k+1} = (2k+1+alpha-x) L_k - (k+alpha) L_{k-1}``. Default ``alpha = 0`` gives the
    ordinary Laguerre polynomials.
    """
    if n < 0:
        raise ValueError("laguerre_l requires n >= 0")
    if n == 0:
        return 1.0
    l_prev = 1.0
    l_cur = 1.0 + alpha - x
    for k in range(1, n):
        l_next = ((2 * k + 1 + alpha - x) * l_cur - (k + alpha) * l_prev) / (k + 1)
        l_prev, l_cur = l_cur, l_next
    return l_cur


def chebyshev_t(n, x):
    """Chebyshev polynomial of the first kind ``T_n(x)`` (weight ``1/sqrt(1-x^2)``).

    ``T_0 = 1``, ``T_1 = x``, ``T_{k+1} = 2x T_k - T_{k-1}``. On ``[-1, 1]``,
    ``T_n(cos theta) = cos(n theta)``, so ``|T_n| <= 1`` there.
    """
    if n < 0:
        raise ValueError("chebyshev_t requires n >= 0")
    if n == 0:
        return 1.0
    t_prev = 1.0
    t_cur = x
    for _ in range(1, n):
        t_next = 2.0 * x * t_cur - t_prev
        t_prev, t_cur = t_cur, t_next
    return t_cur


def chebyshev_u(n, x):
    """Chebyshev polynomial of the second kind ``U_n(x)`` (weight ``sqrt(1-x^2)``).

    ``U_0 = 1``, ``U_1 = 2x``, ``U_{k+1} = 2x U_k - U_{k-1}``. On ``[-1, 1]``,
    ``U_n(cos theta) = sin((n+1) theta) / sin(theta)``.
    """
    if n < 0:
        raise ValueError("chebyshev_u requires n >= 0")
    if n == 0:
        return 1.0
    u_prev = 1.0
    u_cur = 2.0 * x
    for _ in range(1, n):
        u_next = 2.0 * x * u_cur - u_prev
        u_prev, u_cur = u_cur, u_next
    return u_cur
