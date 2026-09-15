import math

from quantforge import (
    laplace_log_evidence,
    bayes_factor,
    posterior_model_probabilities,
)


def close(a, b, tol=1e-4):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_laplace_exact_for_2d_gaussian():
    # Laplace approximation is exact when the log-joint is quadratic
    mu = [1.0, -2.0]
    c = 0.7
    P = [[2.0, 0.5], [0.5, 1.0]]

    def lj(th):
        d0, d1 = th[0] - mu[0], th[1] - mu[1]
        return -0.5 * (P[0][0] * d0 * d0 + 2 * P[0][1] * d0 * d1 + P[1][1] * d1 * d1) + c

    res = laplace_log_evidence(lj, [0.0, 0.0])
    detP = P[0][0] * P[1][1] - P[0][1] * P[1][0]
    analytic = c + math.log(2 * math.pi) - 0.5 * math.log(detP)
    assert close(res["log_evidence"], analytic)
    assert close(res["mode"][0], mu[0]) and close(res["mode"][1], mu[1])


def test_laplace_exact_1d():
    def lj(th):
        return -0.5 * 3.0 * (th[0] - 4.0) ** 2 + 1.2

    r = laplace_log_evidence(lj, [0.0])
    an = 1.2 + 0.5 * math.log(2 * math.pi) - 0.5 * math.log(3.0)
    assert close(r["log_evidence"], an)


def test_bayes_factor():
    assert close(bayes_factor(2.0, 1.0), math.exp(1.0))
    assert close(bayes_factor(1.0, 1.0), 1.0)


def test_posterior_model_probabilities_equal_priors():
    probs = posterior_model_probabilities([math.log(1), math.log(3), math.log(6)])
    assert close(sum(probs), 1.0)
    assert close(probs[0], 0.1) and close(probs[1], 0.3) and close(probs[2], 0.6)


def test_posterior_model_probabilities_with_priors():
    probs = posterior_model_probabilities([0.0, 0.0], priors=[3, 1])
    assert close(probs[0], 0.75) and close(probs[1], 0.25)


def test_model_selection_prefers_prior_near_data():
    data = [4.8, 5.1, 5.2, 4.9, 5.0]

    def make_lj(prior_mu):
        def lj(th):
            mu = th[0]
            ll = sum(-0.5 * (x - mu) ** 2 for x in data)
            lp = -0.5 * (mu - prior_mu) ** 2 / 4.0
            return ll + lp
        return lj

    zA = laplace_log_evidence(make_lj(5.0), [0.0])["log_evidence"]
    zB = laplace_log_evidence(make_lj(0.0), [0.0])["log_evidence"]
    assert zA > zB
    assert bayes_factor(zA, zB) > 1.0
