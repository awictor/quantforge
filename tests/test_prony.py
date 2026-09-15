import cmath
import math

import pytest

from quantforge import prony, prony_reconstruct


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_two_real_exponentials():
    y = [3 * 0.9 ** k + 2 * 0.5 ** k for k in range(20)]
    modes, amps = prony(y, 2)
    mm = sorted(m.real for m in modes)
    assert close(mm[0], 0.5) and close(mm[1], 0.9)
    rec = prony_reconstruct(modes, amps, 20)
    assert max(abs(rec[k] - y[k]) for k in range(20)) < 1e-6


def test_damped_sinusoid():
    y = [math.exp(-0.1 * k) * math.cos(0.5 * k) for k in range(30)]
    modes, amps = prony(y, 2)
    for m in modes:
        assert close(abs(m), math.exp(-0.1), 1e-4)
        assert close(abs(cmath.phase(m)), 0.5, 1e-4)
    rec = prony_reconstruct(modes, amps, 30)
    assert max(abs(rec[k].real - y[k]) for k in range(30)) < 1e-6


def test_single_exponential():
    y = [5 * 0.8 ** k for k in range(10)]
    modes, amps = prony(y, 1)
    assert close(modes[0].real, 0.8) and close(amps[0].real, 5.0)


def test_three_modes():
    y = [2 * 0.95 ** k + 1.5 * 0.7 ** k + 0.5 * 0.4 ** k for k in range(24)]
    modes, _ = prony(y, 3)
    mm = sorted(m.real for m in modes)
    for got, want in zip(mm, (0.4, 0.7, 0.95)):
        assert close(got, want, 1e-4)


def test_undamped_oscillation_on_unit_circle():
    y = [math.cos(0.3 * k) for k in range(20)]
    modes, _ = prony(y, 2)
    for m in modes:
        assert close(abs(m), 1.0, 1e-4)


def test_errors():
    with pytest.raises(ValueError):
        prony([1.0, 2.0, 3.0, 4.0], 0)
    with pytest.raises(ValueError):
        prony([1.0, 2.0], 2)       # need >= 2p samples
