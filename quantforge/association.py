"""Association strength for categorical contingency tables.

The chi-square test says *whether* two categorical variables are associated; these
measures say *how strongly*, on a ``[0, 1]`` scale independent of sample size. Cramer's V
is the general r x c measure; the phi coefficient is its 2x2 special case; Tschuprow's T
and the (bias-corrected) contingency coefficient are alternative normalizations. All are
built from the Pearson chi-square statistic of the table. Pure standard library.
"""

import math


def _chi_square(table):
    """Pearson chi-square statistic and grand total of a contingency table."""
    rows = len(table)
    cols = len(table[0])
    if any(len(r) != cols for r in table):
        raise ValueError("table must be rectangular")
    row_tot = [sum(r) for r in table]
    col_tot = [sum(table[i][j] for i in range(rows)) for j in range(cols)]
    n = sum(row_tot)
    if n == 0:
        raise ValueError("table is empty")
    chi2 = 0.0
    for i in range(rows):
        for j in range(cols):
            expected = row_tot[i] * col_tot[j] / n
            if expected > 0:
                chi2 += (table[i][j] - expected) ** 2 / expected
    return chi2, n, rows, cols


def cramers_v(table):
    """Cramer's V association strength in ``[0, 1]`` for an r x c contingency table.

    ``sqrt(chi2 / (n * min(r-1, c-1)))``. ``0`` for independent variables, ``1`` for a
    perfect association. The general-purpose categorical effect size.
    """
    chi2, n, rows, cols = _chi_square(table)
    k = min(rows - 1, cols - 1)
    if k == 0:
        return 0.0
    return math.sqrt(chi2 / (n * k))


def phi_coefficient(table):
    """Phi coefficient for a 2x2 table: ``sqrt(chi2 / n)`` (equals Cramer's V here).

    Ranges ``[0, 1]`` in magnitude; the signed Pearson-correlation form of a 2x2 table.
    Raises unless the table is 2x2.
    """
    if len(table) != 2 or any(len(r) != 2 for r in table):
        raise ValueError("phi coefficient requires a 2x2 table")
    chi2, n, _, _ = _chi_square(table)
    return math.sqrt(chi2 / n)


def tschuprow_t(table):
    """Tschuprow's T association measure: ``sqrt(chi2 / (n * sqrt((r-1)(c-1))))``.

    Like Cramer's V but reaches ``1`` only for square tables; ``0`` under independence.
    """
    chi2, n, rows, cols = _chi_square(table)
    d = (rows - 1) * (cols - 1)
    if d == 0:
        return 0.0
    return math.sqrt(chi2 / (n * math.sqrt(d)))


def contingency_coefficient(table):
    """Pearson's contingency coefficient ``sqrt(chi2 / (chi2 + n))`` in ``[0, 1)``.

    Always below 1 (it cannot reach it), which is its known limitation; ``0`` under
    independence. Provided for completeness alongside Cramer's V.
    """
    chi2, n, _, _ = _chi_square(table)
    return math.sqrt(chi2 / (chi2 + n))
