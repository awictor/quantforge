import math
import random

from quantforge import extended_kalman_filter
from quantforge.kalman_filter import kalman_filter


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_linear_case_matches_linear_kalman():
    random.seed(0)
    F = [[1.0, 1.0], [0.0, 1.0]]
    H = [[1.0, 0.0]]
    Q = [[1e-3, 0], [0, 1e-3]]
    R = [[0.3]]
    z = [[0.5 * i + random.uniform(-0.4, 0.4)] for i in range(40)]

    def f_lin(x):
        return [x[0] + x[1], x[1]]

    def h_lin(x):
        return [x[0]]

    ek = extended_kalman_filter(z, f_lin, h_lin, Q, R, [0.0, 0.0], [[1.0, 0], [0, 1.0]])
    lk = kalman_filter(z, F, H, Q, R, [0.0, 0.0], [[1.0, 0], [0, 1.0]])
    for i in range(40):
        for j in range(2):
            assert close(ek["filtered_means"][i][j], lk["states"][i][j], 1e-7)


def test_nonlinear_logistic_growth():
    r_g, K = 0.3, 10.0

    def f_log(x):
        return [x[0] + r_g * x[0] * (1 - x[0] / K)]

    def h_id(x):
        return [x[0]]

    true = [0.5]
    for _ in range(39):
        true.append(true[-1] + r_g * true[-1] * (1 - true[-1] / K))
    random.seed(1)
    z = [[t + random.uniform(-0.3, 0.3)] for t in true]
    ek = extended_kalman_filter(z, f_log, h_id, [[1e-4]], [[0.1]], [0.4], [[1.0]])
    est = [m[0] for m in ek["filtered_means"]]
    rmse = math.sqrt(sum((est[i] - true[i]) ** 2 for i in range(40)) / 40)
    assert rmse < 0.2


def test_nonlinear_measurement_model():
    # observe h(x) = x^2, constant state -> recover x
    def f_c(x):
        return [x[0]]

    def h_sq(x):
        return [x[0] * x[0]]

    xt = 3.0
    random.seed(2)
    z = [[xt * xt + random.uniform(-0.5, 0.5)] for _ in range(60)]
    ek = extended_kalman_filter(z, f_c, h_sq, [[1e-5]], [[0.25]], [2.5], [[1.0]])
    assert abs(ek["filtered_means"][-1][0] - xt) < 0.2


def test_covariance_stays_positive():
    def f_c(x):
        return [x[0]]

    def h_sq(x):
        return [x[0] * x[0]]

    random.seed(2)
    z = [[9.0 + random.uniform(-0.5, 0.5)] for _ in range(30)]
    ek = extended_kalman_filter(z, f_c, h_sq, [[1e-5]], [[0.25]], [2.5], [[1.0]])
    for cov in ek["filtered_covariances"]:
        assert cov[0][0] > 0
