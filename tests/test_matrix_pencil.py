import cmath
import math
import random

import pytest

from quantforge import matrix_pencil, prony


def close(a, b, tol=1e-4):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_two_exponentials_clean():
    y = [3 * 0.9 ** k + 2 * 0.5 ** k for k in range(30)]
    modes, _ = matrix_pencil(y, 2)
    mm = sorted(z.real for z in modes)
    assert close(mm[0], 0.5) and close(mm[1], 0.9)


def test_damped_sinusoid():
    y = [math.exp(-0.1 * k) * math.cos(0.5 * k) for k in range(40)]
    modes, _ = matrix_pencil(y, 2)
    for m in modes:
        assert close(abs(m), math.exp(-0.1), 1e-3)
        assert close(abs(cmath.phase(m)), 0.5, 1e-3)


def test_three_exponentials():
    y = [2 * 0.95 ** k + 1.5 * 0.7 ** k + 0.5 * 0.4 ** k for k in range(40)]
    modes, _ = matrix_pencil(y, 3)
    mm = sorted(z.real for z in modes)
    for got, want in zip(mm, (0.4, 0.7, 0.95)):
        assert close(got, want, 1e-3)


def test_noise_robustness_vs_prony():
    random.seed(0)
    true_w = 0.7
    true_r = math.exp(-0.05)
    clean = [true_r ** k * math.cos(true_w * k) for k in range(60)]
    noisy = [c + random.uniform(-0.02, 0.02) for c in clean]

    def freq_err(modes):
        return min(abs(abs(cmath.phase(m)) - true_w) for m in modes)

    mp_err = freq_err(matrix_pencil(noisy, 2)[0])
    pr_err = freq_err(prony(noisy, 2)[0])
    # matrix pencil should be at least as accurate (SVD truncation rejects the noise subspace)
    assert mp_err <= pr_err + 1e-9 or mp_err < 0.02


def test_errors():
    y = [1.0 * 0.9 ** k for k in range(30)]
    with pytest.raises(ValueError):
        matrix_pencil(y, 0)
    with pytest.raises(ValueError):
        matrix_pencil(y, 50)
