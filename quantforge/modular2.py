"""Quadratic residues and multiplicative structure modulo a prime.

The Legendre and Jacobi symbols (is ``a`` a square mod ``p``?), Tonelli-Shanks for the
modular square root, the multiplicative order of an element, and a primitive root
generator. These underpin number-theoretic cryptography and residue arithmetic. Pure
standard library on top of :mod:`quantforge.number_theory` and :mod:`quantforge.modular`.
"""

from .number_theory import is_prime, factorize
from .modular import mod_pow


def legendre_symbol(a, p):
    """Legendre symbol ``(a/p)`` for an odd prime ``p``: 0, 1, or -1.

    ``1`` if ``a`` is a non-zero quadratic residue mod ``p``, ``-1`` if a non-residue,
    ``0`` if ``a`` is divisible by ``p``. Uses Euler's criterion ``a^((p-1)/2) mod p``.
    """
    if p < 3 or not is_prime(p):
        raise ValueError("p must be an odd prime")
    a %= p
    if a == 0:
        return 0
    ls = mod_pow(a, (p - 1) // 2, p)
    return -1 if ls == p - 1 else ls


def jacobi_symbol(a, n):
    """Jacobi symbol ``(a/n)`` for an odd positive ``n``: 0, 1, or -1.

    Generalizes the Legendre symbol to composite (odd) ``n`` by multiplicativity in the
    denominator. Equals the Legendre symbol when ``n`` is prime.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer")
    a %= n
    result = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def tonelli_shanks(a, p):
    """A square root of ``a`` modulo an odd prime ``p`` (Tonelli-Shanks).

    Returns ``r`` with ``r*r == a (mod p)`` (the other root is ``p - r``). Raises if ``a``
    is a non-residue. ``a`` is reduced mod ``p`` first.
    """
    if p < 3 or not is_prime(p):
        raise ValueError("p must be an odd prime")
    a %= p
    if a == 0:
        return 0
    if legendre_symbol(a, p) != 1:
        raise ValueError("a is not a quadratic residue mod p")
    if p % 4 == 3:
        return mod_pow(a, (p + 1) // 4, p)
    # Factor p-1 = q * 2^s with q odd.
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    # Find a non-residue z.
    z = 2
    while legendre_symbol(z, p) != -1:
        z += 1
    m = s
    c = mod_pow(z, q, p)
    t = mod_pow(a, q, p)
    r = mod_pow(a, (q + 1) // 2, p)
    while t != 1:
        # Find the least i, 0 < i < m, with t^(2^i) == 1.
        i = 0
        temp = t
        while temp != 1:
            temp = (temp * temp) % p
            i += 1
        b = mod_pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r


def multiplicative_order(a, n):
    """Multiplicative order of ``a`` modulo ``n``: the least ``k > 0`` with ``a^k == 1``.

    Requires ``gcd(a, n) == 1``. Computes it from the factorization of Euler's totient of
    a prime (here ``n`` is required prime for an exact ``phi = n - 1``); raises otherwise.
    """
    if n < 2 or not is_prime(n):
        raise ValueError("n must be prime")
    a %= n
    if a == 0:
        raise ValueError("a must be coprime to n")
    phi = n - 1
    order = phi
    for prime, _ in factorize(phi):
        while order % prime == 0 and mod_pow(a, order // prime, n) == 1:
            order //= prime
    return order


def primitive_root(p):
    """A primitive root modulo the prime ``p`` (a generator of the multiplicative group).

    Returns the smallest ``g`` whose multiplicative order is ``p - 1``. Every power of a
    primitive root covers all non-zero residues exactly once.
    """
    if p == 2:
        return 1
    if p < 2 or not is_prime(p):
        raise ValueError("p must be prime")
    phi = p - 1
    prime_factors = [q for q, _ in factorize(phi)]
    for g in range(2, p):
        if all(mod_pow(g, phi // q, p) != 1 for q in prime_factors):
            return g
    raise ValueError("no primitive root found")  # unreachable for prime p
