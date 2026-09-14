"""Prime sieves: enumerate, count, and index primes up to a bound.

The Sieve of Eratosthenes lists every prime below a limit in ``O(n log log n)`` -- far
faster than testing each number individually when you need many primes. From it come the
prime-counting function, the n-th prime, and a smallest-prime-factor table for instant
factorization of any number in range. Complements the single-number
:func:`quantforge.is_prime`. Pure standard library.
"""


def primes_up_to(limit):
    """All primes ``<= limit`` by the Sieve of Eratosthenes (a sorted list).

    ``limit < 2`` yields an empty list. Marks multiples starting from each prime's square.
    """
    if limit < 2:
        return []
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    p = 2
    while p * p <= limit:
        if sieve[p]:
            sieve[p * p::p] = bytearray(len(sieve[p * p::p]))
        p += 1
    return [i for i in range(2, limit + 1) if sieve[i]]


def prime_count(limit):
    """Number of primes ``<= limit`` (the prime-counting function ``pi(limit)``)."""
    return len(primes_up_to(limit))


def nth_prime(n):
    """The ``n``-th prime (1-indexed): ``nth_prime(1) == 2``.

    Grows the sieve bound with the prime-number-theorem estimate ``n(ln n + ln ln n)``
    until enough primes are found. ``n >= 1``.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if n < 6:
        return [2, 3, 5, 7, 11][n - 1]
    import math
    # Upper bound (Rosser): p_n < n(ln n + ln ln n) for n >= 6.
    ln = math.log(n)
    limit = int(n * (ln + math.log(ln))) + 3
    primes = primes_up_to(limit)
    while len(primes) < n:
        limit *= 2
        primes = primes_up_to(limit)
    return primes[n - 1]


def smallest_prime_factors(limit):
    """Smallest-prime-factor table for ``0..limit`` (a list; ``spf[k]`` = least prime dividing ``k``).

    ``spf[0] = spf[1] = 0``. With it, any ``k <= limit`` factorizes in ``O(log k)`` by
    repeatedly dividing by ``spf[k]``. Built by a linear-ish sieve.
    """
    if limit < 0:
        raise ValueError("limit must be non-negative")
    spf = list(range(limit + 1))
    if limit >= 0:
        spf[0] = 0
    if limit >= 1:
        spf[1] = 0
    p = 2
    while p * p <= limit:
        if spf[p] == p:                       # p is prime
            for m in range(p * p, limit + 1, p):
                if spf[m] == m:
                    spf[m] = p
        p += 1
    return spf
