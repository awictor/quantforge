import cmath

from quantforge import tf_to_ss, tf_evaluate, tf_dcgain, tf_frequency_response
from quantforge.state_space import lti_impulse_response, dc_gain
from quantforge.eigen_general import eigenvalues_general
from quantforge.polyroots import polynomial_roots


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_first_order_impulse_response():
    A, B, C, D = tf_to_ss([1.0], [1.0, -0.5])   # H = 1/(z - 0.5)
    imp = lti_impulse_response(A, B, C, D, 8)
    assert close(imp[0][0], 0.0)
    for k in range(1, 8):
        assert close(imp[k][0], 0.5 ** (k - 1))


def test_dc_gain_agreement():
    A, B, C, D = tf_to_ss([1.0], [1.0, -0.5])
    assert close(tf_dcgain([1.0], [1.0, -0.5]), dc_gain(A, B, C, D)[0][0])
    assert close(tf_dcgain([1.0], [1.0, -0.5]), 2.0)


def test_poles_are_eigenvalues():
    num, den = [1.0, 0.3], [1.0, -0.9, 0.2]
    A, _, _, _ = tf_to_ss(num, den)
    den_roots = sorted(r.real for r in polynomial_roots(den))
    A_eig = sorted(e.real for e in eigenvalues_general(A))
    for a, b in zip(den_roots, A_eig):
        assert close(a, b)


def test_transfer_function_matches_state_space():
    num, den = [1.0, 0.3], [1.0, -0.9, 0.2]
    A, B, C, D = tf_to_ss(num, den)
    n = len(A)

    def ss_tf(z):
        M = [[complex((z if i == j else 0) - A[i][j]) for j in range(n)] + [complex(B[i][0])]
             for i in range(n)]
        for k in range(n):
            piv = max(range(k, n), key=lambda i: abs(M[i][k]))
            M[k], M[piv] = M[piv], M[k]
            for i in range(k + 1, n):
                f = M[i][k] / M[k][k]
                for j in range(k, n + 1):
                    M[i][j] -= f * M[k][j]
        x = [0j] * n
        for i in range(n - 1, -1, -1):
            x[i] = (M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))) / M[i][i]
        return sum(C[0][j] * x[j] for j in range(n)) + D[0][0]

    for z in (1.5 + 0j, 0.3 + 0.4j, 2.0 + 0j):
        assert abs(tf_evaluate(num, den, z) - ss_tf(z)) < 1e-6


def test_frequency_response():
    fr = tf_frequency_response([1.0], [1.0, -0.5], [0.0, 1.0])
    assert close(fr[0].real, 2.0)
    assert abs(fr[0].imag) < 1e-12
