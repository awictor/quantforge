"""Tests for Bluestein arbitrary-length DFT, cross-checked against a direct DFT and radix-2 FFT."""

import cmath
import random

from quantforge.bluestein import dft, idft
from quantforge.fft import fft


def _direct_dft(x):
    n = len(x)
    return [sum(x[j] * cmath.exp(-2j * cmath.pi * k * j / n) for j in range(n)) for k in range(n)]


def _close(a, b, tol=1e-7):
    return len(a) == len(b) and all(abs(u - v) < tol for u, v in zip(a, b))


def test_fuzz_complex_vs_direct():
    rng = random.Random(131)
    for _ in range(500):
        n = rng.randint(1, 40)
        x = [complex(rng.uniform(-5, 5), rng.uniform(-5, 5)) for _ in range(n)]
        assert _close(dft(x), _direct_dft(x))


def test_fuzz_roundtrip():
    rng = random.Random(132)
    for _ in range(500):
        n = rng.randint(1, 40)
        x = [complex(rng.uniform(-5, 5), rng.uniform(-5, 5)) for _ in range(n)]
        assert _close(idft(dft(x)), x)


def test_fuzz_real_input():
    rng = random.Random(133)
    for _ in range(300):
        n = rng.randint(1, 50)
        x = [rng.uniform(-3, 3) for _ in range(n)]
        assert _close(dft(x), _direct_dft(x))
        assert _close(idft(dft(x)), [complex(v) for v in x])


def test_matches_radix2_on_powers_of_two():
    rng = random.Random(134)
    for p in (1, 2, 4, 8, 16, 32):
        x = [rng.uniform(-2, 2) for _ in range(p)]
        assert _close(dft(x), fft(x))


def test_prime_lengths():
    rng = random.Random(135)
    for n in (3, 5, 7, 11, 13, 17, 19):
        x = [rng.uniform(-1, 1) for _ in range(n)]
        assert _close(dft(x), _direct_dft(x))
        assert _close(idft(dft(x)), [complex(v) for v in x])


def test_dft_of_ones():
    d = dft([1, 1, 1, 1])
    assert abs(d[0] - 4) < 1e-9
    assert all(abs(v) < 1e-9 for v in d[1:])


def test_dft_of_impulse_is_flat():
    d = dft([1, 0, 0])
    assert _close(d, [1, 1, 1])


def test_empty():
    assert dft([]) == []
    assert idft([]) == []


def test_single_element():
    assert _close(dft([5]), [5])
    assert _close(idft([5]), [5])


def test_linearity():
    rng = random.Random(136)
    n = 7  # non-power-of-two
    x = [rng.uniform(-1, 1) for _ in range(n)]
    y = [rng.uniform(-1, 1) for _ in range(n)]
    dx = dft(x)
    dy = dft(y)
    dxy = dft([x[i] + y[i] for i in range(n)])
    assert _close(dxy, [dx[i] + dy[i] for i in range(n)])
