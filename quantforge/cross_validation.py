"""Cross-validation splitters and scoring.

Model-agnostic resampling for out-of-sample evaluation:

  * ``k_fold_indices`` -- partition ``n`` observations into ``k`` disjoint folds
    (optionally shuffled), yielding ``(train_idx, test_idx)`` pairs whose test sets
    exactly tile the data.
  * ``train_test_split`` -- a single deterministic split by fraction.
  * ``cross_val_score`` -- run a caller-supplied fit/predict/score over the folds
    and return the per-fold scores.

A deterministic linear-congruential shuffle keeps splits reproducible per seed.
Pure standard library.
"""


def _lcg_shuffle(n, seed):
    """Fisher-Yates shuffle of range(n) using a deterministic LCG."""
    idx = list(range(n))
    state = seed & 0x7FFFFFFF
    for i in range(n - 1, 0, -1):
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        j = state % (i + 1)
        idx[i], idx[j] = idx[j], idx[i]
    return idx


def k_fold_indices(n, k=5, shuffle=False, seed=1234567):
    """Yield ``(train_indices, test_indices)`` for ``k``-fold cross-validation.

    The ``n`` observations are split into ``k`` contiguous folds (of sizes
    differing by at most one); with ``shuffle=True`` the order is permuted first.
    The test folds are disjoint and cover every index exactly once.
    """
    if k < 2:
        raise ValueError("k must be >= 2")
    if n < k:
        raise ValueError("need at least k observations")
    order = _lcg_shuffle(n, seed) if shuffle else list(range(n))
    # Fold sizes: the first (n mod k) folds get one extra.
    base, extra = divmod(n, k)
    folds = []
    start = 0
    for f in range(k):
        size = base + (1 if f < extra else 0)
        folds.append(order[start:start + size])
        start += size
    for f in range(k):
        test = folds[f]
        train = [i for g in range(k) if g != f for i in folds[g]]
        yield train, test


def train_test_split(n, test_fraction=0.2, shuffle=False, seed=1234567):
    """Single train/test split of ``n`` indices by ``test_fraction``.

    Returns ``(train_indices, test_indices)``. The test set gets
    ``round(n * test_fraction)`` indices (at least 1, at most ``n-1``).
    """
    if not (0.0 < test_fraction < 1.0):
        raise ValueError("test_fraction must be in (0, 1)")
    if n < 2:
        raise ValueError("need at least 2 observations")
    order = _lcg_shuffle(n, seed) if shuffle else list(range(n))
    n_test = min(max(round(n * test_fraction), 1), n - 1)
    test = order[:n_test]
    train = order[n_test:]
    return train, test


def cross_val_score(X, y, fit_fn, score_fn, k=5, shuffle=False, seed=1234567):
    """K-fold cross-validated scores for a caller-supplied model.

    ``fit_fn(X_train, y_train)`` returns a fitted model or predictor; ``score_fn``
    is called as ``score_fn(model, X_test, y_test)`` and returns a scalar. Returns
    the list of ``k`` per-fold scores.
    """
    n = len(y)
    if len(X) != n:
        raise ValueError("X and y must have the same length")
    scores = []
    for train_idx, test_idx in k_fold_indices(n, k, shuffle, seed):
        X_tr = [X[i] for i in train_idx]
        y_tr = [y[i] for i in train_idx]
        X_te = [X[i] for i in test_idx]
        y_te = [y[i] for i in test_idx]
        model = fit_fn(X_tr, y_tr)
        scores.append(score_fn(model, X_te, y_te))
    return scores
