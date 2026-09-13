"""Principal components regression (PCR).

When regressors are many or collinear, OLS variance explodes. PCR first rotates the
(centered) predictors onto their principal components, keeps the top ``k`` (the
directions of greatest variance), regresses ``y`` on just those orthogonal scores, and
maps the fit back to the original predictor space. Dropping low-variance components
regularizes the fit -- trading a little bias for much lower variance -- and with all
``k`` components retained PCR reproduces ordinary least squares exactly. Pure standard
library; builds on the library's PCA.
"""

from .pca import pca


def principal_components_regression(X, y, n_components=None):
    """Principal components regression of ``y`` on the columns of ``X``.

    Centers ``X`` and ``y``, runs PCA on the predictor covariance, keeps the top
    ``n_components`` (default: all), regresses on the component scores, and maps the
    coefficients back to the original variables. Returns a dict with ``coefficients``
    (per original predictor), ``intercept``, ``n_components`` and the
    ``explained_variance`` (cumulative fraction retained). Retaining all components
    reproduces OLS.
    """
    n = len(X)
    if n != len(y):
        raise ValueError("X and y must have equal length")
    if n == 0:
        raise ValueError("need at least one observation")
    p = len(X[0])
    k = p if n_components is None else n_components
    if not (1 <= k <= p):
        raise ValueError("n_components must be in [1, n_features]")

    # Center predictors and response.
    means = [sum(X[t][j] for t in range(n)) / n for j in range(p)]
    ybar = sum(y) / n
    Xc = [[X[t][j] - means[j] for j in range(p)] for t in range(n)]
    yc = [y[t] - ybar for t in range(n)]

    # Predictor covariance (population; scale cancels in the regression).
    cov = [[sum(Xc[t][i] * Xc[t][j] for t in range(n)) / n for j in range(p)]
           for i in range(p)]
    res = pca(cov)
    loadings = res["loadings"]           # loadings[i] = i-th eigenvector (length p)

    # Scores on the top-k components: z_t[i] = Xc[t] . loading_i.
    beta_pc = []
    for i in range(k):
        vec = loadings[i]
        z = [sum(Xc[t][j] * vec[j] for j in range(p)) for t in range(n)]
        zz = sum(zt * zt for zt in z)
        if zz <= 0:
            beta_pc.append(0.0)
        else:
            beta_pc.append(sum(z[t] * yc[t] for t in range(n)) / zz)

    # Map back: beta = sum_i beta_pc[i] * loading_i.
    beta = [0.0] * p
    for i in range(k):
        vec = loadings[i]
        b = beta_pc[i]
        for j in range(p):
            beta[j] += b * vec[j]

    intercept = ybar - sum(beta[j] * means[j] for j in range(p))
    explained = res["cumulative_explained"][k - 1]
    return {"coefficients": beta, "intercept": intercept,
            "n_components": k, "explained_variance": explained}
