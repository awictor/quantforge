import math

from quantforge import remez
from quantforge.remez import _poly_eval
from quantforge.chebyshev import chebyshev_fit, chebyshev_eval


def _maxerr(f, coeffs, a, b, N=5000):
    return max(abs(f(a + (b - a) * i / N) - _poly_eval(coeffs, a + (b - a) * i / N))
               for i in range(N + 1))


def test_linear_is_exact():
    r = remez(lambda x: 2 * x + 1, -1.0, 1.0, 2)
    assert r["error"] < 1e-9


def test_exp_error_decreases_with_degree():
    prev = 1e9
    for deg in (1, 2, 3, 4):
        r = remez(math.exp, -1.0, 1.0, deg)
        me = _maxerr(math.exp, r["coeffs"], -1.0, 1.0)
        assert abs(me - r["error"]) < 1e-3        # reported error matches actual max error
        assert r["error"] < prev
        prev = r["error"]
    # degree-3 minimax error for exp on [-1,1] is ~0.0055 (known)
    assert abs(remez(math.exp, -1.0, 1.0, 3)["error"] - 0.00553) < 5e-4


def test_minimax_beats_least_squares_in_max_error():
    cases = [(math.exp, -1.0, 1.0, 4), (math.sin, 0.0, math.pi, 5),
             (lambda x: 1.0 / (1.0 + x), 0.0, 2.0, 4)]
    for f, a, b, deg in cases:
        rm = remez(f, a, b, deg)
        cf = chebyshev_fit(f, a, b, deg)
        cheb_max = max(abs(f(a + (b - a) * i / 2000) - chebyshev_eval(cf, a, b, a + (b - a) * i / 2000))
                       for i in range(2001))
        assert rm["error"] <= cheb_max + 1e-9


def test_equioscillation():
    deg = 4
    r = remez(math.exp, -1.0, 1.0, deg)
    N = 5000
    errs = [math.exp(-1 + 2 * i / N) - _poly_eval(r["coeffs"], -1 + 2 * i / N) for i in range(N + 1)]
    sign_changes = sum(1 for i in range(1, len(errs)) if errs[i - 1] * errs[i] < 0)
    assert sign_changes >= deg + 1        # n+1 roots -> n+2 alternating extrema
