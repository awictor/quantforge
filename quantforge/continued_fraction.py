"""Continued-fraction expansion and best rational approximation of a real number.

Every real number has a continued-fraction expansion ``a0 + 1/(a1 + 1/(a2 + ...))``
whose truncations -- the *convergents* -- are the best rational approximations of their
size: no fraction with a smaller denominator is closer. This is the tool behind rational
approximation of constants, gear-ratio and calendar design, and recovering a simple
fraction from a noisy decimal. ``cf_expansion`` computes the coefficients,
``convergents`` the successive fractions, and ``best_rational`` the closest fraction
whose denominator stays under a bound. Pure standard library.
"""

import math


def cf_expansion(x, max_terms=64, tol=1e-15):
    """Continued-fraction coefficients ``[a0, a1, a2, ...]`` of a real ``x``.

    Repeatedly takes the integer part and inverts the remainder. Stops after
    ``max_terms`` terms or when the fractional remainder falls below ``tol`` (a rational
    ``x`` terminates). ``a0`` may be negative or zero; every later term is a positive
    integer.
    """
    terms = []
    a0 = math.floor(x)
    terms.append(int(a0))
    frac = x - a0
    for _ in range(max_terms - 1):
        if frac <= tol:
            break
        inv = 1.0 / frac
        a = math.floor(inv)
        terms.append(int(a))
        frac = inv - a
    return terms


def convergents(terms):
    """Successive convergents ``[(p0, q0), (p1, q1), ...]`` from CF coefficients.

    Uses the standard recurrence ``p_k = a_k p_{k-1} + p_{k-2}``,
    ``q_k = a_k q_{k-1} + q_{k-2}``. Each ``(p, q)`` is the best rational approximation
    with denominator ``<= q`` -- the fractions converge to the value, alternating above
    and below it.
    """
    if not terms:
        raise ValueError("need at least one term")
    p_prev, p_prev2 = 1, 0        # p_{-1}=1, p_{-2}=0
    q_prev, q_prev2 = 0, 1        # q_{-1}=0, q_{-2}=1
    out = []
    for a in terms:
        p = a * p_prev + p_prev2
        q = a * q_prev + q_prev2
        out.append((p, q))
        p_prev2, p_prev = p_prev, p
        q_prev2, q_prev = q_prev, q
    return out


def best_rational(x, max_denominator=1000000):
    """Closest fraction ``(p, q)`` to ``x`` with ``q <= max_denominator``.

    Returns the convergent (or the appropriate *semiconvergent*) with the largest
    admissible denominator -- the best rational approximation under the bound, matching
    the classic Stern-Brocot / ``limit_denominator`` result. ``q >= 1`` always.
    """
    if max_denominator < 1:
        raise ValueError("max_denominator must be >= 1")
    terms = cf_expansion(x, max_terms=128)
    conv = convergents(terms)
    best = (conv[0][0], conv[0][1])
    for k in range(len(conv)):
        p, q = conv[k]
        if q > max_denominator:
            # Try a semiconvergent using the previous convergent to get closer within bound.
            if k >= 1:
                p_prev, q_prev = conv[k - 1]
                # largest t with q_prev + t*q_{prev-1-step}... use the CF term reduction:
                # semiconvergent denominator q_{k-2} + t*q_{k-1} <= max_denominator.
                p_prev2, q_prev2 = (conv[k - 2] if k >= 2 else (terms[0] and 1 or 1, 0))
                if q_prev != 0:
                    t = (max_denominator - q_prev2) // q_prev
                    if t > 0:
                        cand_p = p_prev2 + t * p_prev
                        cand_q = q_prev2 + t * q_prev
                        if cand_q <= max_denominator and cand_q >= 1:
                            # accept the semiconvergent only if it beats the last convergent
                            if abs(cand_p / cand_q - x) < abs(best[0] / best[1] - x):
                                best = (cand_p, cand_q)
            break
        best = (p, q)
    return best
