"""Prime sieves: enumerate, count, index, factor."""

import random

import pytest

from quantforge import primes_up_to, prime_count, nth_prime, smallest_prime_factors
from quantforge import is_prime, factorize


def test_sieve_vs_is_prime():
    for lim in (0, 1, 2, 3, 50, 1000):
        assert set(primes_up_to(lim)) == {n for n in range(2, lim + 1) if is_prime(n)}


def test_prime_count():
    assert prime_count(10) == 4
    assert prime_count(100) == 25
    assert prime_count(1000) == 168


def test_nth_prime():
    for n, p in [(1, 2), (6, 13), (25, 97), (100, 541), (1000, 7919), (10000, 104729)]:
        assert nth_prime(n) == p


def test_smallest_prime_factors():
    spf = smallest_prime_factors(10000)
    rng = random.Random(1)
    for _ in range(2000):
        k = rng.randint(2, 10000)
        fac = {}
        m = k
        while m > 1:
            p = spf[m]
            fac[p] = fac.get(p, 0) + 1
            m //= p
        assert sorted(fac.items()) == factorize(k)


def test_validation():
    with pytest.raises(ValueError):
        nth_prime(0)
    with pytest.raises(ValueError):
        smallest_prime_factors(-1)
