from quantforge import (
    solve_discrete_lyapunov,
    solve_continuous_lyapunov,
    controllability_gramian,
)


def _mm(X, Y):
    return [[sum(X[i][k] * Y[k][j] for k in range(len(Y))) for j in range(len(Y[0]))]
            for i in range(len(X))]


def _tr(X):
    return [[X[i][j] for i in range(len(X))] for j in range(len(X[0]))]


def close(a, b, tol=1e-8):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_discrete_residual_and_symmetry():
    A = [[0.5, 0.1], [0.0, 0.4]]
    Q = [[1.0, 0.0], [0.0, 1.0]]
    P = solve_discrete_lyapunov(A, Q)
    APAt = _mm(_mm(A, P), _tr(A))
    res = [[APAt[i][j] - P[i][j] + Q[i][j] for j in range(2)] for i in range(2)]
    assert max(abs(res[i][j]) for i in range(2) for j in range(2)) < 1e-8
    assert close(P[0][1], P[1][0])


def test_discrete_matches_series():
    A = [[0.5, 0.1], [0.0, 0.4]]
    Q = [[1.0, 0.0], [0.0, 1.0]]
    P = solve_discrete_lyapunov(A, Q)
    Pit = [[0.0, 0.0], [0.0, 0.0]]
    Ak = [[1.0, 0.0], [0.0, 1.0]]
    for _ in range(200):
        term = _mm(_mm(Ak, Q), _tr(Ak))
        Pit = [[Pit[i][j] + term[i][j] for j in range(2)] for i in range(2)]
        Ak = _mm(A, Ak)
    for i in range(2):
        for j in range(2):
            assert close(P[i][j], Pit[i][j], 1e-6)


def test_continuous_residual():
    A = [[-1.0, 0.2], [0.0, -2.0]]
    Q = [[1.0, 0.0], [0.0, 1.0]]
    P = solve_continuous_lyapunov(A, Q)
    lhs = _mm(A, P)
    lhs = [[lhs[i][j] + _mm(P, _tr(A))[i][j] + Q[i][j] for j in range(2)] for i in range(2)]
    assert max(abs(lhs[i][j]) for i in range(2) for j in range(2)) < 1e-8


def test_scalar_closed_forms():
    assert close(solve_discrete_lyapunov([[0.6]], [[2.0]])[0][0], 2.0 / (1 - 0.36))
    assert close(solve_continuous_lyapunov([[-1.5]], [[3.0]])[0][0], -3.0 / (2 * -1.5))


def test_controllability_gramian_positive_definite():
    A = [[0.5, 0.1], [0.0, 0.4]]
    B = [[0.0], [1.0]]
    W = controllability_gramian(A, B, discrete=True)
    assert W[0][0] > 0
    assert W[0][0] * W[1][1] - W[0][1] * W[1][0] > 0
