"""Probabilistic-forecast calibration: Brier decomposition, reliability, ECE.

For binary-outcome probability forecasts, Murphy's decomposition splits the Brier
score into three interpretable pieces,

    Brier = reliability - resolution + uncertainty,

where *reliability* measures how far the observed frequency in each forecast group
drifts from the forecast itself (0 = perfectly calibrated), *resolution* rewards
forecasts that separate outcomes away from the base rate, and *uncertainty* is the
irreducible variance of the outcome ``obar (1 - obar)``. Grouping forecasts by
identical value makes the identity exact; binning continuous forecasts makes it
approximate. Also provides the reliability curve (calibration diagram) and the
expected calibration error. Pure standard library.
"""


def _bin_forecasts(forecasts, outcomes, n_bins):
    """Group into bins; return per-group (mean_forecast, obs_freq, count)."""
    n = len(forecasts)
    if n != len(outcomes):
        raise ValueError("forecasts and outcomes must have equal length")
    if n == 0:
        raise ValueError("need at least one forecast")
    if any(not (0.0 <= p <= 1.0) for p in forecasts):
        raise ValueError("forecasts must lie in [0, 1]")
    if any(o not in (0, 1, 0.0, 1.0) for o in outcomes):
        raise ValueError("outcomes must be 0 or 1")

    if n_bins is None:
        # Exact grouping by identical forecast value.
        groups = {}
        for p, o in zip(forecasts, outcomes):
            groups.setdefault(p, []).append(o)
        keys = sorted(groups)
        return [(k, sum(groups[k]) / len(groups[k]), len(groups[k])) for k in keys]

    # Equal-width bins on [0, 1]; edge 1.0 falls in the last bin.
    buckets = [[] for _ in range(n_bins)]
    fsum = [0.0] * n_bins
    for p, o in zip(forecasts, outcomes):
        idx = int(p * n_bins)
        if idx == n_bins:
            idx = n_bins - 1
        buckets[idx].append(o)
        fsum[idx] += p
    out = []
    for k in range(n_bins):
        cnt = len(buckets[k])
        if cnt == 0:
            continue
        out.append((fsum[k] / cnt, sum(buckets[k]) / cnt, cnt))
    return out


def brier_decomposition(forecasts, outcomes, n_bins=None):
    """Murphy decomposition of the Brier score into reliability/resolution/uncertainty.

    Returns a dict with ``reliability``, ``resolution``, ``uncertainty``, ``brier``
    (the reconstructed ``reliability - resolution + uncertainty``), and ``base_rate``.
    With ``n_bins=None`` forecasts are grouped by identical value and the identity is
    exact; with an integer ``n_bins`` they are put into equal-width bins on ``[0, 1]``
    and the reconstruction is approximate.
    """
    groups = _bin_forecasts(forecasts, outcomes, n_bins)
    n = len(forecasts)
    base = sum(outcomes) / n
    reliability = sum(cnt * (f - o) ** 2 for f, o, cnt in groups) / n
    resolution = sum(cnt * (o - base) ** 2 for f, o, cnt in groups) / n
    uncertainty = base * (1.0 - base)
    return {
        "reliability": reliability,
        "resolution": resolution,
        "uncertainty": uncertainty,
        "brier": reliability - resolution + uncertainty,
        "base_rate": base,
    }


def reliability_curve(forecasts, outcomes, n_bins=10):
    """Calibration diagram: per-bin mean forecast vs observed frequency.

    Returns a list of ``(mean_forecast, observed_frequency, count)`` tuples, one per
    non-empty equal-width bin. A perfectly calibrated forecaster lies on the diagonal
    ``observed == forecast``.
    """
    return _bin_forecasts(forecasts, outcomes, n_bins)


def expected_calibration_error(forecasts, outcomes, n_bins=10):
    """Expected calibration error: count-weighted mean ``|forecast - observed|`` per bin.

    Zero for a perfectly calibrated forecaster; the standard scalar summary of a
    reliability diagram's departure from the diagonal.
    """
    groups = _bin_forecasts(forecasts, outcomes, n_bins)
    n = len(forecasts)
    return sum(cnt * abs(f - o) for f, o, cnt in groups) / n
