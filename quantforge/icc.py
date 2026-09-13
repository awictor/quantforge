"""Intraclass correlation coefficient (ICC) for rater / measurement reliability.

The ICC measures how strongly measurements on the same subject resemble each other --
what fraction of the total variance is *between* subjects rather than measurement
noise. Given a subjects x raters table it decomposes the variance with a two-way ANOVA
and forms the Shrout-Fleiss coefficients:

  * ``ICC(1)``   -- one-way: each subject rated by different raters (rater is noise).
  * ``ICC(2,1)`` / ``ICC(2,k)`` -- two-way random: raters are a random sample; single
    rating vs the mean of ``k`` raters (absolute agreement).
  * ``ICC(3,1)`` / ``ICC(3,k)`` -- two-way mixed: these raters are the whole population
    (consistency, rater bias not counted as error).

``1`` = perfect agreement, ``0`` = no better than chance, and it can go slightly
negative. Pure standard library.
"""


def _anova_components(data):
    n = len(data)           # subjects
    k = len(data[0])        # raters
    grand = sum(sum(row) for row in data) / (n * k)

    row_means = [sum(row) / k for row in data]
    col_means = [sum(data[i][j] for i in range(n)) / n for j in range(k)]

    # Sum of squares.
    ss_total = sum((data[i][j] - grand) ** 2 for i in range(n) for j in range(k))
    ss_rows = k * sum((rm - grand) ** 2 for rm in row_means)          # between subjects
    ss_cols = n * sum((cm - grand) ** 2 for cm in col_means)          # between raters
    ss_err = ss_total - ss_rows - ss_cols                            # residual

    ms_rows = ss_rows / (n - 1)
    ms_cols = ss_cols / (k - 1)
    ms_err = ss_err / ((n - 1) * (k - 1))
    ms_within = (ss_cols + ss_err) / (n * (k - 1))                    # one-way error
    return n, k, ms_rows, ms_cols, ms_err, ms_within


def icc(data):
    """Shrout-Fleiss intraclass correlation coefficients.

    ``data`` is a list of rows, one per subject, each a length-``k`` list of the
    ``k`` raters' scores (a balanced subjects x raters table). Returns a dict with
    ``icc1``, ``icc2_1``, ``icc2_k``, ``icc3_1`` and ``icc3_k``. Requires at least two
    subjects and two raters.
    """
    rows = [list(r) for r in data]
    n = len(rows)
    if n < 2:
        raise ValueError("need at least 2 subjects")
    k = len(rows[0])
    if k < 2:
        raise ValueError("need at least 2 raters")
    if any(len(r) != k for r in rows):
        raise ValueError("all subjects must have the same number of raters")

    n, k, msr, msc, mse, msw = _anova_components(rows)

    icc1 = (msr - msw) / (msr + (k - 1) * msw) if (msr + (k - 1) * msw) != 0 else 0.0
    icc1k = (msr - msw) / msr if msr != 0 else 0.0

    denom_21 = msr + (k - 1) * mse + k * (msc - mse) / n
    icc21 = (msr - mse) / denom_21 if denom_21 != 0 else 0.0
    denom_2k = msr + (msc - mse) / n
    icc2k = (msr - mse) / denom_2k if denom_2k != 0 else 0.0

    icc31 = (msr - mse) / (msr + (k - 1) * mse) if (msr + (k - 1) * mse) != 0 else 0.0
    icc3k = (msr - mse) / msr if msr != 0 else 0.0

    return {
        "icc1": icc1,
        "icc1k": icc1k,
        "icc2_1": icc21,
        "icc2_k": icc2k,
        "icc3_1": icc31,
        "icc3_k": icc3k,
    }
