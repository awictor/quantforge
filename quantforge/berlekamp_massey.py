"""Berlekamp-Massey: the shortest linear recurrence a sequence satisfies (mod a prime).

Given the first terms of a sequence over a prime field, Berlekamp-Massey finds the shortest
linear recurrence ``s_i = sum_j c_j s_{i-1-j} (mod p)`` that reproduces them -- the same
object a linear-feedback shift register realizes. It runs in ``O(n^2)`` and, given at least
``2L`` terms of a sequence with a length-``L`` recurrence, recovers it exactly. Useful for
guessing a closed form, decoding, or extrapolating a modular sequence. Pure standard
library.
"""


def berlekamp_massey(sequence, mod):
    """Return the coefficients of the shortest recurrence for ``sequence`` modulo ``mod``.

    ``mod`` must be prime. The returned list ``c`` has length ``L`` (the recurrence order)
    and satisfies ``sequence[i] == sum_j c[j] * sequence[i-1-j] (mod mod)`` for all
    ``i >= L``. An all-zero sequence returns ``[]`` (order 0).
    """
    s = [x % mod for x in sequence]
    n = len(s)
    ls, cur = [], []          # last (before shortening) and current recurrence coeff lists
    lf = 0                    # index at which ls was captured
    ld = 0                    # discrepancy at that capture
    for i in range(n):
        # discrepancy d = s[i] - (recurrence's prediction of s[i])
        t = 0
        for j in range(len(cur)):
            t = (t + cur[j] * s[i - 1 - j]) % mod
        d = (s[i] - t) % mod
        if d == 0:
            continue
        if not cur:
            cur = [0] * (i + 1)
            lf = i
            ld = d
            continue
        k = (d * pow(ld, mod - 2, mod)) % mod
        # c = cur + k * x^(i-lf) applied to (unit - ls) at the right offset
        c = [0] * (i - lf - 1) + [k] + [(-k * x) % mod for x in ls]
        if len(c) < len(cur):
            c += [0] * (len(cur) - len(c))
        for j in range(len(cur)):
            c[j] = (c[j] + cur[j]) % mod
        if i - len(cur) > lf - len(ls):
            ls, lf, ld = cur[:], i, d
        cur = c
    return [x % mod for x in cur]


def _next_term(coeffs, history, mod):
    return sum(coeffs[j] * history[-1 - j] for j in range(len(coeffs))) % mod


def berlekamp_massey_next(sequence, mod, count=1):
    """Extend ``sequence`` by ``count`` terms using its shortest recurrence (mod ``mod``).

    Runs Berlekamp-Massey, then rolls the recurrence forward. Returns the list of the next
    ``count`` terms. Raises if the sequence is too short to determine any recurrence and a
    prediction is requested.
    """
    if count < 0:
        raise ValueError("count must be non-negative")
    coeffs = berlekamp_massey(sequence, mod)
    hist = [x % mod for x in sequence]
    out = []
    if not coeffs:
        # zero recurrence: the sequence continues as all zeros
        for _ in range(count):
            out.append(0)
            hist.append(0)
        return out
    for _ in range(count):
        nxt = _next_term(coeffs, hist, mod)
        out.append(nxt)
        hist.append(nxt)
    return out
