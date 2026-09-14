"""Integer number theory: primality, factorization, gcd/lcm, totient.

The building blocks of exact integer work: a deterministic Miller-Rabin primality test
(exact for all 64-bit integers), Pollard's rho for fast factorization of large
composites, and the derived quantities -- the prime factorization, gcd/lcm, and Euler's
totient. These underpin cryptographic sizing, modular arithmetic, and any exact-integer
computation. Pure standard library (Python big integers are exact, so there is no
overflow).
"""


def gcd(a, b):
    """Greatest common divisor of ``a`` and ``b`` (Euclid's algorithm, non-negative)."""
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def lcm(a, b):
    """Least common multiple of ``a`` and ``b`` (``0`` if either is zero)."""
    if a == 0 or b == 0:
        return 0
    return abs(a // gcd(a, b) * b)


def is_prime(n):
    """Deterministic Miller-Rabin primality test, exact for all ``n < 3.3 * 10^24``.

    Uses a fixed set of witness bases that is proven to give no false positives below
    that bound (well past 64-bit), so the answer is exact, not probabilistic, for any
    integer arising in normal use. ``O(k log^3 n)``.
    """
    if n < 2:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small_primes:
        if n % p == 0:
            return n == p
    # Write n-1 = d * 2^r.
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in small_primes:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _pollard_rho(n):
    """One nontrivial factor of a composite ``n`` by Pollard's rho (Brent's variant)."""
    if n % 2 == 0:
        return 2
    # Deterministic pseudo-random walk; vary the constant if a cycle fails to split.
    for c in range(1, n):
        f = lambda x: (x * x + c) % n
        x = y = 2
        d = 1
        while d == 1:
            x = f(x)
            y = f(f(y))
            d = gcd(abs(x - y), n)
        if d != n:
            return d
    return n   # unreachable for composite n


def factorize(n):
    """Prime factorization of ``n`` as a sorted list of ``(prime, exponent)`` pairs.

    Combines trial division by small primes with Pollard's rho for large factors and
    Miller-Rabin to certify primality, so it factors numbers far beyond what naive trial
    division reaches. ``n`` must be ``>= 1`` (``1`` has an empty factorization).
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if n == 1:
        return []
    factors = {}

    def _add(p):
        factors[p] = factors.get(p, 0) + 1

    stack = [n]
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            _add(m)
            continue
        d = _pollard_rho(m)
        stack.append(d)
        stack.append(m // d)
    return sorted(factors.items())


def divisors(n):
    """All positive divisors of ``n`` in sorted order (from its factorization)."""
    if n < 1:
        raise ValueError("n must be >= 1")
    divs = [1]
    for p, e in factorize(n):
        divs = [d * p ** k for d in divs for k in range(e + 1)]
    return sorted(divs)


def euler_totient(n):
    """Euler's totient ``phi(n)``: the count of integers in ``[1, n]`` coprime to ``n``.

    Computed from the factorization as ``n * prod (1 - 1/p)`` over distinct primes ``p``.
    ``phi(1) = 1``.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    result = n
    for p, _ in factorize(n):
        result -= result // p
    return result
