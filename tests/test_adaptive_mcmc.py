import math

from quantforge import adaptive_metropolis, slice_sample
from quantforge.mcmc import sample_mean, sample_cov


def _corr_logp():
    mean = [1.0, -1.0]
    S = [[2.0, 0.8], [0.8, 1.0]]
    det = S[0][0] * S[1][1] - S[0][1] * S[1][0]
    P = [[S[1][1] / det, -S[0][1] / det], [-S[1][0] / det, S[0][0] / det]]

    def lp(x):
        dx0, dx1 = x[0] - mean[0], x[1] - mean[1]
        return -0.5 * (P[0][0] * dx0 * dx0 + 2 * P[0][1] * dx0 * dx1 + P[1][1] * dx1 * dx1)

    return lp


def test_adaptive_metropolis_correlated_gaussian():
    lp = _corr_logp()
    res = adaptive_metropolis(lp, [0.0, 0.0], 30000, seed=1, burn_in=5000)
    m = sample_mean(res["samples"])
    c = sample_cov(res["samples"])
    assert abs(m[0] - 1) < 0.15 and abs(m[1] + 1) < 0.15
    assert abs(c[0][0] - 2) < 0.4 and abs(c[1][1] - 1) < 0.3 and abs(c[0][1] - 0.8) < 0.3
    assert 0.1 < res["accept_rate"] < 0.7


def test_slice_sampler_1d_gaussian():
    mu, sig = 2.0, 1.5
    s = slice_sample(lambda x: -0.5 * ((x - mu) / sig) ** 2, 0.0, 20000, w=3.0, seed=2, burn_in=1000)
    sm = sum(s) / len(s)
    sv = sum((v - sm) ** 2 for v in s) / len(s)
    assert abs(sm - mu) < 0.1
    assert abs(sv - sig * sig) < 0.3


def test_slice_sampler_visits_both_modes():
    def lp(x):
        a = math.exp(-0.5 * ((x + 3) / 0.7) ** 2)
        b = math.exp(-0.5 * ((x - 3) / 0.7) ** 2)
        return math.log(a + b + 1e-300)

    s = slice_sample(lp, 0.0, 20000, w=2.0, seed=3, burn_in=2000)
    neg = sum(1 for v in s if v < 0)
    assert abs(sum(s) / len(s)) < 0.5
    assert neg > 2000 and len(s) - neg > 2000


def test_reproducible():
    lp = _corr_logp()
    a = adaptive_metropolis(lp, [0.0, 0.0], 2000, seed=9, burn_in=100)
    b = adaptive_metropolis(lp, [0.0, 0.0], 2000, seed=9, burn_in=100)
    assert a["samples"] == b["samples"]
    lp1 = lambda x: -0.5 * x ** 2
    assert slice_sample(lp1, 0.0, 1000, seed=9) == slice_sample(lp1, 0.0, 1000, seed=9)
