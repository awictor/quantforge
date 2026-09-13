"""Recursive least squares: online OLS updated one observation at a time.

Batch OLS re-solves ``(X'X)^{-1} X'y`` from scratch. Recursive least squares (RLS)
updates the coefficient vector and the inverse covariance ``P = (X'X)^{-1}`` in
``O(k^2)`` per new observation using the Sherman-Morrison identity -- no re-inversion,
no storing the full history. With ``lambda = 1`` it converges to the exact batch OLS
fit; a forgetting factor ``lambda < 1`` down-weights old data, tracking slowly-varying
coefficients (adaptive filtering, time-varying betas). Pure standard library.
"""


class RecursiveLeastSquares:
    """Online multivariate least squares with an optional forgetting factor.

    Feed observations one at a time with :meth:`update`; read the current fit from
    :attr:`beta`. ``n_features`` is the regressor count (include a constant 1 in each
    ``x`` for an intercept). ``forgetting`` in ``(0, 1]`` down-weights past data
    (``1`` = ordinary growing-window OLS); ``delta`` sets the prior ``P = delta * I``
    (large = diffuse prior). With ``forgetting = 1`` the estimate matches batch OLS
    once enough points have arrived.
    """

    def __init__(self, n_features, forgetting=1.0, delta=1e6):
        if n_features < 1:
            raise ValueError("n_features must be >= 1")
        if not (0.0 < forgetting <= 1.0):
            raise ValueError("forgetting must be in (0, 1]")
        self.k = n_features
        self.lam = forgetting
        self.beta = [0.0] * n_features
        # P = delta * I (inverse of the regularized X'X).
        self.P = [[delta if i == j else 0.0 for j in range(n_features)]
                  for i in range(n_features)]
        self.n = 0

    def update(self, x, y):
        """Incorporate one observation ``(x, y)`` and return the updated coefficients.

        ``x`` is a length-``n_features`` regressor row, ``y`` the scalar response.
        """
        if len(x) != self.k:
            raise ValueError("x has wrong length")
        k = self.k
        lam = self.lam
        P = self.P

        # Px = P @ x.
        Px = [sum(P[i][j] * x[j] for j in range(k)) for i in range(k)]
        # Denominator lambda + x' P x.
        denom = lam + sum(x[i] * Px[i] for i in range(k))
        # Gain g = Px / denom.
        g = [Px[i] / denom for i in range(k)]
        # Prediction error.
        pred = sum(self.beta[i] * x[i] for i in range(k))
        err = y - pred
        # Update beta.
        self.beta = [self.beta[i] + g[i] * err for i in range(k)]
        # Update P = (P - g (Px)') / lambda.
        self.P = [[(P[i][j] - g[i] * Px[j]) / lam for j in range(k)] for i in range(k)]
        self.n += 1
        return list(self.beta)

    def predict(self, x):
        """Predict the response for a regressor row ``x`` from the current fit."""
        if len(x) != self.k:
            raise ValueError("x has wrong length")
        return sum(self.beta[i] * x[i] for i in range(self.k))


def recursive_least_squares(X, y, forgetting=1.0, delta=1e6):
    """Fit RLS over a whole dataset and return the final coefficient vector.

    Convenience wrapper: streams the rows of ``X`` (each already including any
    intercept column) through :class:`RecursiveLeastSquares`. With ``forgetting = 1``
    the result matches batch OLS on the same design.
    """
    n = len(y)
    if n != len(X):
        raise ValueError("X and y must have equal length")
    if n == 0:
        raise ValueError("need at least one observation")
    rls = RecursiveLeastSquares(len(X[0]), forgetting=forgetting, delta=delta)
    for t in range(n):
        rls.update(list(map(float, X[t])), float(y[t]))
    return rls.beta
