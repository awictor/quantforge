import random

from quantforge import lti_simulate, lti_step_response, lti_impulse_response, dc_gain


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_scalar_first_order_step():
    a = 0.5
    A, B, C, D = [[a]], [[1.0]], [[1.0]], [[0.0]]
    ys = lti_step_response(A, B, C, D, 10)
    for k in range(10):
        assert close(ys[k][0], (1 - a ** k) / (1 - a))


def test_step_converges_to_dc_gain():
    A, B, C, D = [[0.5]], [[1.0]], [[1.0]], [[0.0]]
    g = dc_gain(A, B, C, D)[0][0]
    assert close(g, 2.0)
    assert close(lti_step_response(A, B, C, D, 40)[-1][0], g, 1e-6)


def test_impulse_markov_parameters():
    a = 0.5
    A, B, C, D = [[a]], [[1.0]], [[1.0]], [[0.0]]
    imp = lti_impulse_response(A, B, C, D, 8)
    assert close(imp[0][0], 0.0)                 # h_0 = D = 0
    for k in range(1, 8):
        assert close(imp[k][0], a ** (k - 1))    # h_k = C A^{k-1} B


def test_direct_feedthrough():
    ys = lti_step_response([[0.3]], [[1.0]], [[1.0]], [[2.0]], 5)
    assert close(ys[0][0], 2.0)                  # y_0 = C*0 + D*1


def test_linearity():
    random.seed(0)
    A, B, C, D = [[0.5, 0.1], [0.0, 0.4]], [[0.0], [1.0]], [[1.0, 0.0]], [[0.0]]
    u = [[random.gauss(0, 1)] for _ in range(20)]
    u2 = [[2 * x[0]] for x in u]
    y1 = lti_simulate(A, B, C, D, u)
    y2 = lti_simulate(A, B, C, D, u2)
    for k in range(20):
        assert close(y2[k][0], 2 * y1[k][0])


def test_impulse_sum_is_step():
    A, B, C, D = [[0.5, 0.1], [0.0, 0.4]], [[0.0], [1.0]], [[1.0, 0.0]], [[0.0]]
    imp = [r[0] for r in lti_impulse_response(A, B, C, D, 50)]
    step = [r[0] for r in lti_step_response(A, B, C, D, 50)]
    for k in range(50):
        assert close(sum(imp[:k + 1]), step[k])
