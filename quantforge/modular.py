"""Modular arithmetic: extended GCD, inverse, CRT, discrete logarithm.

The core operations of modular (clock) arithmetic beyond ``pow(a, b, m)``: the extended
Euclidean algorithm and the modular inverse it yields, the Chinese Remainder Theorem for
stitching congruences with coprime moduli into one, and a baby-step/giant-step discrete
logarithm. These underpin RSA-style key arithmetic, hashing, and any exact computation
in a residue ring. Pure standard library on exact Python integers.
"""


def extended_gcd(a, b):
    """Extended Euclidean algorithm: return ``(g, x, y)`` with ``a*x + b*y == g = gcd(a, b)``.

    The Bezout coefficients ``x, y`` are the basis of the modular inverse and the CRT.
    """
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t


def mod_inverse(a, m):
    """Modular inverse ``a^{-1} mod m``: the ``x`` with ``a*x == 1 (mod m)``.

    Exists iff ``gcd(a, m) == 1``; raises ``ValueError`` otherwise. Result is in
    ``[0, m)``.
    """
    if m <= 0:
        raise ValueError("modulus must be positive")
    g, x, _ = extended_gcd(a % m, m)
    if g != 1:
        raise ValueError("no inverse: %d and %d are not coprime" % (a, m))
    return x % m


def chinese_remainder(remainders, moduli):
    """Chinese Remainder Theorem: the unique ``x`` in ``[0, M)`` matching all congruences.

    Given ``x == remainders[i] (mod moduli[i])`` with pairwise-coprime ``moduli``,
    returns ``(x, M)`` where ``M`` is the product of the moduli. Raises if the moduli are
    not pairwise coprime.
    """
    if len(remainders) != len(moduli):
        raise ValueError("remainders and moduli must have equal length")
    if not moduli:
        raise ValueError("need at least one congruence")
    M = 1
    for m in moduli:
        if m <= 0:
            raise ValueError("moduli must be positive")
        M *= m
    x = 0
    for r, m in zip(remainders, moduli):
        Mi = M // m
        g, inv, _ = extended_gcd(Mi % m, m)
        if g != 1:
            raise ValueError("moduli are not pairwise coprime")
        x += r * Mi * (inv % m)
    return x % M, M


def mod_pow(base, exp, mod):
    """Modular exponentiation ``base^exp mod mod`` (supports negative ``exp`` via inverse).

    A thin, explicit wrapper over fast binary exponentiation; a negative exponent inverts
    the base first (requires ``gcd(base, mod) == 1``).
    """
    if mod <= 0:
        raise ValueError("modulus must be positive")
    if exp < 0:
        base = mod_inverse(base, mod)
        exp = -exp
    return pow(base, exp, mod)


def discrete_log(base, target, mod):
    """Smallest non-negative ``x`` with ``base^x == target (mod mod)``, or ``None``.

    Baby-step giant-step: ``O(sqrt(mod))`` time and space. Searches exponents in
    ``[0, mod)``. ``mod`` must be positive; ``target`` is reduced mod ``mod``.
    """
    if mod <= 0:
        raise ValueError("modulus must be positive")
    target %= mod
    if mod == 1:
        return 0
    n = 1
    while n * n < mod:
        n += 1                       # n = ceil(sqrt(mod))
    # Baby steps: base^j for j in [0, n).
    table = {}
    e = 1
    for j in range(n):
        table.setdefault(e, j)       # keep the smallest j for each value
        e = (e * base) % mod
    # Giant steps: target * (base^-n)^i.
    factor = mod_inverse(pow(base, n, mod), mod)
    gamma = target
    for i in range(n + 1):
        if gamma in table:
            x = i * n + table[gamma]
            if pow(base, x, mod) == target:
                return x
        gamma = (gamma * factor) % mod
    return None
