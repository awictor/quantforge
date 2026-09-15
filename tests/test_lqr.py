from quantforge import lqr, dare, controllability_matrix, is_controllable
from quantforge.eigen_general import eigenvalues_general


def close(a, b, tol=1e-5):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _mm(X, Y):
    return [[sum(X[i][k] * Y[k][j] for k in range(len(Y))) for j in range(len(Y[0]))]
            for i in range(len(X))]


def _tr(X):
    return [[X[i][j] for i in range(len(X))] for j in range(len(X[0]))]


def test_scalar_lqr_closed_form():
    a, b, q, r = 1.1, 1.0, 1.0, 1.0
    res = lqr([[a]], [[b]], [[q]], [[r]])
    p = res["P"][0][0]
    K = res["K"][0][0]
    S = r + p * b * b
    assert abs(a * a * p - a * a * p * p * b * b / S + q - p) < 1e-8   # DARE residual
    assert close(K, b * p * a / S)
    assert abs(a - b * K) < 1.0                                        # closed-loop stable


def test_double_integrator_stabilizes():
    A = [[1.0, 1.0], [0.0, 1.0]]
    B = [[0.0], [1.0]]
    Q = [[1.0, 0.0], [0.0, 1.0]]
    R = [[1.0]]
    res = lqr(A, B, Q, R)
    K = res["K"]
    BK = [[B[i][0] * K[0][j] for j in range(2)] for i in range(2)]
    Acl = [[A[i][j] - BK[i][j] for j in range(2)] for i in range(2)]
    for e in eigenvalues_general(Acl):
        assert abs(complex(e)) < 1.0


def test_dare_fixed_point():
    A = [[1.0, 1.0], [0.0, 1.0]]
    B = [[0.0], [1.0]]
    Q = [[1.0, 0.0], [0.0, 1.0]]
    R = [[1.0]]
    P = dare(A, B, Q, R)
    At, Bt = _tr(A), _tr(B)
    BtP = _mm(Bt, P)
    S = R[0][0] + _mm(BtP, B)[0][0]
    AtP = _mm(At, P)
    term = _mm(_mm(_mm(AtP, B), [[1 / S]]), _mm(BtP, A))
    rhs = [[_mm(AtP, A)[i][j] - term[i][j] + Q[i][j] for j in range(2)] for i in range(2)]
    assert max(abs(rhs[i][j] - P[i][j]) for i in range(2) for j in range(2)) < 1e-6


def test_controllability():
    A = [[1.0, 1.0], [0.0, 1.0]]
    B = [[0.0], [1.0]]
    assert is_controllable(A, B)
    # diagonal system where B touches only one mode is uncontrollable
    assert not is_controllable([[2.0, 0.0], [0.0, 3.0]], [[1.0], [0.0]])


def test_controllability_matrix_shape():
    A = [[1.0, 1.0], [0.0, 1.0]]
    B = [[0.0], [1.0]]
    C = controllability_matrix(A, B)
    assert len(C) == 2 and len(C[0]) == 2   # n x (n*m) = 2 x 2
