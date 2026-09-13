"""Adaptive Gauss-Kronrod quadrature (G7-K15) with an embedded error estimate.

The 15-point Kronrod rule reuses the 7-point Gauss rule's nodes and adds eight more,
so a single set of function evaluations yields both a high-order estimate (K15) and,
by comparison with the embedded G7 estimate, a cheap local error bound. Subdividing
the interval where that error is largest (a global adaptive scheme) concentrates work
on the hard parts of the integrand. This is the workhorse behind QUADPACK's ``QAG``.
Pure standard library.
"""

# Kronrod-15 abscissae (nonnegative half) and weights, and the Gauss-7 weights,
# on the reference interval [-1, 1] (QUADPACK QK15 constants).
_XGK = [
    0.991455371120813, 0.949107912342759, 0.864864423359769,
    0.741531185599394, 0.586087235467691, 0.405845151377397,
    0.207784955007898, 0.000000000000000,
]
_WGK = [
    0.022935322010529, 0.063092092629979, 0.104790010322250,
    0.140653259715525, 0.169004726639267, 0.190350578064785,
    0.204432940075298, 0.209482141084728,
]
# Gauss-7 weights (apply to the odd-indexed Kronrod abscissae 1,3,5,7).
_WG = [
    0.129484966168870, 0.279705391489277, 0.381830050505119,
    0.417959183673469,
]


def _gk15(f, a, b):
    """One G7-K15 panel over [a, b]: returns (kronrod_estimate, error_estimate)."""
    center = 0.5 * (a + b)
    half = 0.5 * (b - a)
    fc = f(center)
    # Kronrod sum and Gauss sum.
    resk = _WGK[7] * fc
    resg = _WG[3] * fc
    fv = [0.0] * 7
    for j in range(7):
        x = half * _XGK[j]
        f1 = f(center - x)
        f2 = f(center + x)
        fsum = f1 + f2
        resk += _WGK[j] * fsum
        if j % 2 == 1:                       # odd Kronrod index = Gauss node
            resg += _WG[j // 2] * fsum
    resk *= half
    resg *= half
    return resk, abs(resk - resg)


def gauss_kronrod(f, a, b, tol=1e-10, max_intervals=1000):
    """Adaptive Gauss-Kronrod (G7-K15) integral of ``f`` over ``[a, b]``.

    Starts with one panel and repeatedly bisects the panel with the largest local
    error estimate until the total estimated error falls below ``tol`` (or
    ``max_intervals`` panels are used). Returns the integral. Concentrates
    evaluations where the integrand is hardest, so it handles peaks and mild
    endpoint behaviour that a fixed rule would miss.
    """
    r, e = _gk15(f, a, b)
    panels = [(a, b, r, e)]
    total = r
    total_err = e
    while total_err > tol and len(panels) < max_intervals:
        # Split the worst panel.
        wi = max(range(len(panels)), key=lambda i: panels[i][3])
        pa, pb, pr, pe = panels.pop(wi)
        mid = 0.5 * (pa + pb)
        r1, e1 = _gk15(f, pa, mid)
        r2, e2 = _gk15(f, mid, pb)
        panels.append((pa, mid, r1, e1))
        panels.append((mid, pb, r2, e2))
        total = total - pr + r1 + r2
        total_err = total_err - pe + e1 + e2
    return total
