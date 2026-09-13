"""McNemar and Cochran's Q tests for paired / repeated binary outcomes.

Paired-binary analogues of the paired t-test and Friedman test:

  * ``mcnemar_test`` -- two paired binary measurements (e.g. a diagnostic before and
    after, or two classifiers on the same items). Only the *discordant* pairs carry
    information: ``b`` = (0 then 1), ``c`` = (1 then 0). The test asks whether
    ``b == c``; it returns the exact two-sided binomial p-value (``b ~ Binomial(b+c,
    1/2)``) plus the continuity-corrected chi-square approximation.
  * ``cochran_q_test`` -- ``k`` binary treatments measured on the same ``b`` blocks
    (the binary special case of Friedman). Tests equal success rates across treatments
    with a chi-square-``(k-1)`` statistic.

Pure standard library.
"""

from .distributions import binomial_cdf, chi2_cdf


def mcnemar_test(table=None, b=None, c=None):
    """McNemar's test for a paired 2x2 table.

    Provide either the 2x2 ``table`` ``[[a, b], [c, d]]`` (a/d concordant, b/c
    discordant) or the two discordant counts ``b`` and ``c`` directly. Returns a dict
    with the discordant counts, the exact two-sided binomial ``p_value``, and the
    continuity-corrected chi-square statistic ``chi2_cc`` with its ``p_value_chi2``.
    """
    if table is not None:
        if len(table) != 2 or any(len(r) != 2 for r in table):
            raise ValueError("table must be 2x2")
        b = table[0][1]
        c = table[1][0]
    if b is None or c is None:
        raise ValueError("provide a 2x2 table or both b and c")
    if b < 0 or c < 0:
        raise ValueError("counts must be non-negative")

    n = b + c
    if n == 0:
        return {"b": b, "c": c, "p_value": 1.0, "chi2_cc": 0.0, "p_value_chi2": 1.0}

    # Exact two-sided binomial on min(b, c) ~ Binomial(n, 1/2).
    k = min(b, c)
    p_exact = min(2.0 * binomial_cdf(k, n, 0.5), 1.0)

    # Continuity-corrected chi-square (Edwards).
    chi2_cc = (abs(b - c) - 1.0) ** 2 / n if n > 0 else 0.0
    if abs(b - c) < 1.0:
        chi2_cc = 0.0
    p_chi2 = 1.0 - chi2_cdf(chi2_cc, 1)
    return {"b": b, "c": c, "p_value": p_exact,
            "chi2_cc": chi2_cc, "p_value_chi2": p_chi2}


def cochran_q_test(blocks):
    """Cochran's Q test for ``k`` binary treatments over ``b`` blocks.

    ``blocks`` is a sequence of rows, each a length-``k`` sequence of 0/1 outcomes for
    one block across the treatments. Returns a dict with the ``statistic`` Q, ``df``
    (``k - 1``) and the chi-square ``p_value``. Blocks whose outcomes are all-0 or
    all-1 contribute nothing (as in the standard formulation).
    """
    rows = [list(r) for r in blocks]
    b = len(rows)
    if b < 2:
        raise ValueError("need at least 2 blocks")
    k = len(rows[0])
    if k < 2:
        raise ValueError("need at least 2 treatments")
    if any(len(r) != k for r in rows):
        raise ValueError("all blocks must have the same number of treatments")
    if any(v not in (0, 1, 0.0, 1.0) for r in rows for v in r):
        raise ValueError("outcomes must be 0 or 1")

    col_totals = [sum(rows[i][j] for i in range(b)) for j in range(k)]
    row_totals = [sum(r) for r in rows]
    grand = sum(col_totals)

    sum_col_sq = sum(cj * cj for cj in col_totals)
    sum_row = sum(row_totals)
    sum_row_sq = sum(ri * ri for ri in row_totals)

    denom = k * sum_row - sum_row_sq
    if denom == 0:
        q = 0.0
    else:
        q = (k - 1) * (k * sum_col_sq - grand * grand) / denom
    df = k - 1
    p = 1.0 - chi2_cdf(q, df)
    return {"statistic": q, "df": df, "p_value": p}
