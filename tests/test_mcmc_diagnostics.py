import pytest

from quantforge import (
    gelman_rubin,
    integrated_autocorrelation_time,
    effective_sample_size,
)
from quantforge.pcg import PCG32
from quantforge.particle_filter import pcg_gaussian


def close(a, b, tol=0.15):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _ar1(phi, n, seed=7):
    r = PCG32(seed)
    x = [0.0]
    for _ in range(n):
        x.append(phi * x[-1] + pcg_gaussian(r, 0, 1))
    return x[1000:]


def test_white_noise_tau_and_ess():
    r = PCG32(1)
    wn = [pcg_gaussian(r, 0, 1) for _ in range(20000)]
    assert abs(integrated_autocorrelation_time(wn) - 1.0) < 0.3
    assert effective_sample_size(wn) > 15000


def test_ar1_autocorrelation_time():
    for phi in (0.5, 0.8):
        tau = integrated_autocorrelation_time(_ar1(phi, 40000))
        expect = (1 + phi) / (1 - phi)
        assert close(tau, expect, 0.2)


def test_ar1_effective_sample_size():
    x = _ar1(0.8, 40000)
    expect_tau = (1 + 0.8) / (1 - 0.8)
    assert close(effective_sample_size(x), len(x) / expect_tau, 0.25)


def _chain(seed, off=0.0, n=5000):
    r = PCG32(seed)
    return [pcg_gaussian(r, 0, 1) + off for _ in range(n)]


def test_gelman_rubin_convergent():
    chains = [_chain(s) for s in (1, 2, 3, 4)]
    assert abs(gelman_rubin(chains) - 1.0) < 0.05


def test_gelman_rubin_divergent():
    div = [_chain(s, off) for s, off in ((1, 0), (2, 5), (3, 10), (4, 15))]
    assert gelman_rubin(div) > 3.0


def test_ess_never_exceeds_n():
    assert effective_sample_size([1.0, 2.0, 1.0, 2.0, 1.0, 2.0]) <= 6


def test_errors():
    with pytest.raises(ValueError):
        gelman_rubin([[1, 2, 3]])          # single chain
    with pytest.raises(ValueError):
        integrated_autocorrelation_time([1.0])   # too short
