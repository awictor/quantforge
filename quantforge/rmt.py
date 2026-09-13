"""Random-matrix-theory denoising of a correlation matrix (eigenvalue clipping).

A sample correlation matrix estimated from ``T`` observations of ``N`` assets has
most of its eigenvalues buried in random noise. Under the null of independent
returns those noise eigenvalues fall below the Marchenko-Pastur upper edge

    lambda_plus = (1 + sqrt(N / T))^2

(for unit-variance data). Eigenvalues above the edge carry signal. Clipping replaces
the sub-edge (noise) eigenvalues with their common average -- preserving the trace,
so the matrix stays a correlation matrix -- while keeping the signal eigenvalues
intact. The denoised matrix is far better-conditioned for portfolio optimization.
Pure standard library on top of the Jacobi eigensolver.
"""

from .pca import jacobi_eigen


def marchenko_pastur_edge(n_assets, n_obs):
    """Upper edge ``(1 + sqrt(N/T))^2`` of the Marchenko-Pastur spectrum.

    Eigenvalues of a noise correlation matrix (unit variances) lie below this; those
    above it are candidate signal. Requires ``n_obs >= n_assets`` for a full-rank
    sample.
    """
    if n_assets < 1 or n_obs < 1:
        raise ValueError("n_assets and n_obs must be positive")
    q = n_assets / n_obs
    return (1.0 + q ** 0.5) ** 2


def clip_correlation_eigenvalues(correlation, n_obs):
    """Denoise a correlation matrix by clipping sub-Marchenko-Pastur eigenvalues.

    Eigen-decomposes ``correlation``, replaces every eigenvalue below the
    Marchenko-Pastur edge with the average of those noise eigenvalues (keeping the
    signal eigenvalues), and rebuilds the matrix, then rescales the diagonal back to
    exactly one. The trace is preserved and the result is a valid, better-conditioned
    correlation matrix. ``n_obs`` is the number of observations used to estimate
    ``correlation``.
    """
    n = len(correlation)
    if n == 0 or any(len(row) != n for row in correlation):
        raise ValueError("correlation must be a non-empty square matrix")
    if n_obs < 1:
        raise ValueError("n_obs must be positive")
    eigenvalues, eigenvectors = jacobi_eigen(correlation)
    edge = marchenko_pastur_edge(n, n_obs)

    noise = [e for e in eigenvalues if e < edge]
    if noise:
        avg_noise = sum(noise) / len(noise)
        clipped = [e if e >= edge else avg_noise for e in eigenvalues]
    else:
        clipped = list(eigenvalues)

    # Rebuild C = sum_k lambda_k v_k v_k'. eigenvectors[k] pairs with eigenvalues[k].
    out = [[0.0] * n for _ in range(n)]
    for k in range(n):
        v = eigenvectors[k]
        lam = clipped[k]
        for i in range(n):
            lv = lam * v[i]
            for j in range(n):
                out[i][j] += lv * v[j]
    # Rescale to unit diagonal (numerical cleanup; clipping preserves the trace).
    d = [out[i][i] ** 0.5 if out[i][i] > 0 else 1.0 for i in range(n)]
    for i in range(n):
        for j in range(n):
            out[i][j] /= (d[i] * d[j])
    return out
