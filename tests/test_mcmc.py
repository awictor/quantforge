import math

from quantforge import (
    metropolis_hastings,
    hamiltonian_monte_carlo,
    sample_mean,
    sample_cov,
)


def test_metropolis_recovers_1d_gaussian():
    mu, sig = 2.0, 1.5
    res = metropolis_hastings(lambda x: -0.5 * ((x[0] - mu) / sig) ** 2,
                              [0.0], 20000, step=2.0, seed=1, burn_in=2000)
    m = sample_mean(res["samples"])[0]
    v = sample_cov(res["samples"])[0][0]
    assert abs(m - mu) < 0.1
    assert abs(v - sig * sig) < 0.4
    assert 0.1 < res["accept_rate"] < 0.9


def test_hmc_recovers_1d_gaussian():
    mu, sig = 2.0, 1.5
    res = hamiltonian_monte_carlo(lambda x: -0.5 * ((x[0] - mu) / sig) ** 2,
                                  [0.0], 5000, step=0.3, n_leapfrog=15, seed=2, burn_in=500)
    m = sample_mean(res["samples"])[0]
    v = sample_cov(res["samples"])[0][0]
    assert abs(m - mu) < 0.1
    assert abs(v - sig * sig) < 0.4
    assert res["accept_rate"] > 0.6


def _corr_gaussian_logp():
    mean = [1.0, -1.0]
    S = [[2.0, 0.8], [0.8, 1.0]]
    det = S[0][0] * S[1][1] - S[0][1] * S[1][0]
    P = [[S[1][1] / det, -S[0][1] / det], [-S[1][0] / det, S[0][0] / det]]

    def logp(x):
        dx0 = x[0] - mean[0]
        dx1 = x[1] - mean[1]
        q = P[0][0] * dx0 * dx0 + 2 * P[0][1] * dx0 * dx1 + P[1][1] * dx1 * dx1
        return q * (-0.5)

    return logp, mean, S


def test_metropolis_2d_correlated():
    logp, mean, S = _corr_gaussian_logp()
    res = metropolis_hastings(logp, [0.0, 0.0], 40000, step=1.5, seed=3, burn_in=4000)
    m = sample_mean(res["samples"])
    c = sample_cov(res["samples"])
    assert abs(m[0] - 1) < 0.15 and abs(m[1] + 1) < 0.15
    assert abs(c[0][0] - 2) < 0.4 and abs(c[1][1] - 1) < 0.3 and abs(c[0][1] - 0.8) < 0.3


def test_hmc_2d_correlated():
    logp, mean, S = _corr_gaussian_logp()
    res = hamiltonian_monte_carlo(logp, [0.0, 0.0], 6000, step=0.25, n_leapfrog=20,
                                  seed=4, burn_in=1000)
    m = sample_mean(res["samples"])
    c = sample_cov(res["samples"])
    assert abs(m[0] - 1) < 0.2 and abs(m[1] + 1) < 0.2
    assert abs(c[0][1] - 0.8) < 0.35


def test_reproducible_with_seed():
    logp = lambda x: -0.5 * x[0] ** 2
    a = metropolis_hastings(logp, [0.0], 1000, seed=9)
    b = metropolis_hastings(logp, [0.0], 1000, seed=9)
    assert a["samples"] == b["samples"]
