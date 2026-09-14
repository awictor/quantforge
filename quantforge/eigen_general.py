"""General (non-symmetric) eigenvalues via the Faddeev-LeVerrier method.

Where ``jacobi_eigen`` handles symmetric matrices (real eigenvalues, orthogonal
eigenvectors), a general real matrix can have complex eigenvalues in conjugate pairs. The
Faddeev-LeVerrier recurrence builds the characteristic polynomial's coefficients directly
from traces of matrix powers -- and as a byproduct yields the determinant and, when
non-singular, the inverse. Rooting that polynomial (Durand-Kerner) gives every eigenvalue,
real or complex. Pure standard library.
"""

from .polyroots import polynomial_roots


def _matmul(A, B):
    n = len(A)
    p = len(B[0])
    m = len(B)
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]


def _trace(A):
    return sum(A[i][i] for i in range(len(A)))


def characteristic_polynomial(A):
    """Coefficients of ``det(xI - A)`` (monic, highest degree first) by Faddeev-LeVerrier.

    Returns ``[1, c_{n-1}, ..., c_0]`` -- the characteristic polynomial with leading
    coefficient 1. Also the input to :func:`polynomial_roots` for the eigenvalues. ``A``
    must be square.
    """
    n = len(A)
    if n == 0 or any(len(row) != n for row in A):
        raise ValueError("matrix must be square and non-empty")
    # Faddeev-LeVerrier: M_1 = I, c_1 = -tr(A M_1); M_{k+1} = A M_k + c_k I.
    coeffs = [1.0]
    I = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    M = [row[:] for row in I]
    for k in range(1, n + 1):
        AM = _matmul(A, M)
        c = -_trace(AM) / k
        coeffs.append(c)
        # M_{k+1} = A M_k + c I
        M = [[AM[i][j] + (c if i == j else 0.0) for j in range(n)] for i in range(n)]
    return coeffs


def eigenvalues_general(A, tol=1e-10):
    """All eigenvalues of a general real (or complex) square matrix.

    Builds the characteristic polynomial via Faddeev-LeVerrier and finds its roots with
    Durand-Kerner, so complex-conjugate eigenvalue pairs are returned correctly. Real
    eigenvalues come back as Python ``float``; genuinely complex ones as ``complex``.
    Sorted by real part then imaginary part.
    """
    coeffs = characteristic_polynomial(A)
    roots = polynomial_roots(coeffs)
    out = []
    for r in roots:
        if abs(r.imag) < tol:
            out.append(r.real)
        else:
            out.append(complex(r.real, r.imag))
    out.sort(key=lambda z: (z.real, z.imag) if isinstance(z, complex) else (z, 0.0))
    return out


def determinant_from_charpoly(A):
    """Determinant of ``A`` from its characteristic polynomial's constant term.

    ``det(A) = (-1)^n * c_0`` where ``c_0`` is the constant term of ``det(xI - A)``. A
    cheap exact-by-construction cross-check against the LU determinant.
    """
    n = len(A)
    coeffs = characteristic_polynomial(A)
    return (-1) ** n * coeffs[-1]
