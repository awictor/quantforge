"""Integer number theory: primality, factorization, gcd/lcm, totient."""

import random

import pytest

from quantforge import gcd, lcm, is_prime, factorize, divisors, euler_totient


def _brute_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def test_primality_matches_brute_force():
    assert all(is_prime(n) == _brute_prime(n) for n in range(2001))


def test_primality_hard_cases():
    assert is_prime(2 ** 61 - 1)          # Mersenne prime
    assert not is_prime(561)              # Carmichael number
    assert not is_prime(1)
    assert is_prime(10 ** 18 + 9)


def test_factorize_reconstructs_and_certifies():
    for n in [2, 60, 997 * 991, 2 ** 10 * 3 ** 5 * 7, 999999999989]:
        f = factorize(n)
        prod = 1
        for p, e in f:
            prod *= p ** e
        assert prod == n
        assert all(is_prime(p) for p, _ in f)
    assert factorize(60) == [(2, 2), (3, 1), (5, 1)]
    assert factorize(1) == []


def test_factorize_large_semiprime():
    p, q = 1000000007, 1000000009
    assert factorize(p * q) == [(p, 1), (q, 1)]


def test_totient_matches_brute_force():
    def brute(n):
        return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)
    assert all(euler_totient(n) == brute(n) for n in range(1, 500))


def test_divisors():
    assert divisors(28) == [1, 2, 4, 7, 14, 28]
    assert sum(divisors(28)[:-1]) == 28          # 28 is perfect
    assert divisors(1) == [1]


def test_gcd_lcm_identity():
    rng = random.Random(1)
    for _ in range(1000):
        a = rng.randint(1, 10 ** 6)
        b = rng.randint(1, 10 ** 6)
        assert gcd(a, b) * lcm(a, b) == a * b
    assert lcm(0, 5) == 0


def test_validation():
    with pytest.raises(ValueError):
        factorize(0)
    with pytest.raises(ValueError):
        euler_totient(0)
    with pytest.raises(ValueError):
        divisors(-3)
