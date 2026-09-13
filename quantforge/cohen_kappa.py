"""Cohen's kappa, weighted kappa, and Fleiss' kappa for categorical agreement.

Where the ICC measures agreement on *continuous* ratings, kappa measures it on
*categorical* ones, correcting the raw agreement for the amount expected by chance:

    kappa = (p_observed - p_expected) / (1 - p_expected).

  * ``cohen_kappa`` -- two raters, nominal categories.
  * ``weighted_kappa`` -- two raters, *ordinal* categories, penalizing disagreements
    by how far apart they are (``linear`` or ``quadratic`` weights) so a near-miss
    counts less than a wild one.
  * ``fleiss_kappa`` -- ``m`` raters per subject (fixed number), from a subjects x
    categories count table.

``1`` = perfect agreement, ``0`` = chance level, negative = worse than chance. Pure
standard library.
"""


def _confusion(rater_a, rater_b):
    if len(rater_a) != len(rater_b):
        raise ValueError("raters must have equal length")
    if len(rater_a) == 0:
        raise ValueError("need at least one item")
    cats = sorted(set(rater_a) | set(rater_b))
    index = {c: i for i, c in enumerate(cats)}
    k = len(cats)
    m = [[0 for _ in range(k)] for _ in range(k)]
    for a, b in zip(rater_a, rater_b):
        m[index[a]][index[b]] += 1
    return m, k


def cohen_kappa(rater_a, rater_b):
    """Cohen's kappa for two raters over nominal categories.

    ``rater_a`` and ``rater_b`` are equal-length label sequences. Returns the
    chance-corrected agreement in ``[-1, 1]``.
    """
    m, k = _confusion(rater_a, rater_b)
    n = sum(sum(row) for row in m)
    p_obs = sum(m[i][i] for i in range(k)) / n
    row_tot = [sum(m[i]) for i in range(k)]
    col_tot = [sum(m[i][j] for i in range(k)) for j in range(k)]
    p_exp = sum(row_tot[i] * col_tot[i] for i in range(k)) / (n * n)
    if p_exp == 1.0:
        return 1.0
    return (p_obs - p_exp) / (1.0 - p_exp)


def weighted_kappa(rater_a, rater_b, weights="linear"):
    """Weighted kappa for two raters over ordinal categories.

    Disagreements are penalized by category distance: ``weights="linear"`` uses
    ``|i - j| / (k - 1)`` and ``"quadratic"`` uses ``(i - j)^2 / (k - 1)^2``. The
    labels must be sortable into their ordinal order. Quadratic weighting is the common
    choice and, for a square table, coincides with an ICC-style measure.
    """
    m, k = _confusion(rater_a, rater_b)
    if k == 1:
        return 1.0
    n = sum(sum(row) for row in m)
    row_tot = [sum(m[i]) for i in range(k)]
    col_tot = [sum(m[i][j] for i in range(k)) for j in range(k)]

    def w(i, j):
        if weights == "linear":
            return abs(i - j) / (k - 1)
        if weights == "quadratic":
            return (i - j) ** 2 / ((k - 1) ** 2)
        raise ValueError("weights must be 'linear' or 'quadratic'")

    obs = sum(w(i, j) * m[i][j] for i in range(k) for j in range(k)) / n
    exp = sum(w(i, j) * row_tot[i] * col_tot[j]
              for i in range(k) for j in range(k)) / (n * n)
    if exp == 0:
        return 1.0
    return 1.0 - obs / exp


def fleiss_kappa(table):
    """Fleiss' kappa for ``m`` raters per subject over fixed categories.

    ``table`` is a list of rows, one per subject, each giving the count of raters
    assigning each category (every row sums to the same ``m``). Returns the
    chance-corrected multi-rater agreement.
    """
    rows = [list(r) for r in table]
    n = len(rows)
    if n < 2:
        raise ValueError("need at least 2 subjects")
    k = len(rows[0])
    if any(len(r) != k for r in rows):
        raise ValueError("all subjects must span the same categories")
    m = sum(rows[0])
    if m < 2:
        raise ValueError("need at least 2 ratings per subject")
    if any(sum(r) != m for r in rows):
        raise ValueError("every subject must have the same number of ratings")

    # Per-subject agreement.
    p_i = [(sum(c * c for c in rows[i]) - m) / (m * (m - 1)) for i in range(n)]
    p_bar = sum(p_i) / n
    # Category marginals.
    p_j = [sum(rows[i][j] for i in range(n)) / (n * m) for j in range(k)]
    p_e = sum(pj * pj for pj in p_j)
    if p_e == 1.0:
        return 1.0
    return (p_bar - p_e) / (1.0 - p_e)
