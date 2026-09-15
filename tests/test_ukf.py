import math
import random

from quantforge import unscented_kalman_filter as ukf
from quantforge.kalman_filter import kalman_filter


def close(a, b, tol=1e-5):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_linear_case_matches_linear_kalman():
    random.seed(0)
    F = [[1.0, 1.0], [0.0, 1.0]]
    H = [[1.0, 0.0]]
    Q = [[1e-3, 0], [0, 1e-3]]
    R = [[0.3]]
    z = [[0.5 * i + random.uniform(-0.4, 0.4)] for i in range(40)]
    u = ukf(z, lambda x: [x[0] + x[1], x[1]], lambda x: [x[0]], Q, R,
            [0.0, 0.0], [[1.0, 0], [0, 1.0]], alpha=1.0, beta=0.0, kappa=0.0)
    lk = kalman_filter(z, F, H, Q, R, [0.0, 0.0], [[1.0, 0], [0, 1.0]])
    for i in range(40):
        for j in range(2):
            assert close(u["filtered_means"][i][j], lk["states"][i][j], 1e-6)


def test_logistic_growth_tracking():
    r_g, K = 0.3, 10.0
    true = [0.5]
    for _ in range(39):
        true.append(true[-1] + r_g * true[-1] * (1 - true[-1] / K))
    random.seed(1)
    z = [[t + random.uniform(-0.3, 0.3)] for t in true]
    u = ukf(z, lambda x: [x[0] + r_g * x[0] * (1 - x[0] / K)], lambda x: [x[0]],
            [[1e-4]], [[0.1]], [0.4], [[1.0]])
    est = [m[0] for m in u["filtered_means"]]
    rmse = math.sqrt(sum((est[i] - true[i]) ** 2 for i in range(40)) / 40)
    assert rmse < 0.2


def test_nonlinear_measurement_recovers_state():
    random.seed(2)
    z = [[9.0 + random.uniform(-0.5, 0.5)] for _ in range(60)]
    u = ukf(z, lambda x: [x[0]], lambda x: [x[0] * x[0]],
            [[1e-5]], [[0.25]], [2.5], [[1.0]])
    assert abs(u["filtered_means"][-1][0] - 3.0) < 0.2


def test_covariance_positive():
    random.seed(2)
    z = [[9.0 + random.uniform(-0.5, 0.5)] for _ in range(30)]
    u = ukf(z, lambda x: [x[0]], lambda x: [x[0] * x[0]],
            [[1e-5]], [[0.25]], [2.5], [[1.0]])
    for cov in u["filtered_covariances"]:
        assert cov[0][0] > 0


def test_agrees_with_ekf_on_mild_nonlinearity():
    from quantforge import extended_kalman_filter as ekf

    random.seed(3)
    z = [[math.sin(0.1 * i) + random.uniform(-0.2, 0.2)] for i in range(50)]
    u = ukf(z, lambda x: [x[0] + 0.1 * x[1], x[1] - 0.05 * math.sin(x[0])],
            lambda x: [x[0]], [[1e-3, 0], [0, 1e-3]], [[0.1]],
            [0.0, 0.0], [[1.0, 0], [0, 1.0]])
    e = ekf(z, lambda x: [x[0] + 0.1 * x[1], x[1] - 0.05 * x[0].sin()],
            lambda x: [x[0]], [[1e-3, 0], [0, 1e-3]], [[0.1]],
            [0.0, 0.0], [[1.0, 0], [0, 1.0]])
    for j in range(2):
        assert abs(u["filtered_means"][-1][j] - e["filtered_means"][-1][j]) < 0.05
