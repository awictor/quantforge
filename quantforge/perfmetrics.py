"""Performance and drawdown statistics for a return series.

Standard track-record measures computed from a sequence of periodic returns
(simple, not log): the Sharpe and Sortino ratios (annualized), maximum
drawdown of the cumulative-return curve, the Calmar ratio, and hit rate /
profit factor. All pure standard library.
"""

import math
from typing import Sequence


def _mean(xs):
    return sum(xs) / len(xs)


def _std(xs, ddof=1):
    n = len(xs)
    if n - ddof <= 0:
        raise ValueError("need more observations than ddof")
    m = _mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (n - ddof))


def cumulative_return(returns: Sequence[float]) -> float:
    """Total compounded return over the series, ``prod(1 + r) - 1``."""
    growth = 1.0
    for r in returns:
        growth *= (1.0 + r)
    return growth - 1.0


def annualized_return(returns: Sequence[float], periods_per_year=252) -> float:
    """Geometric (compound) annualized return.

    ``(prod(1 + r))^{periods_per_year / n} - 1`` -- the constant per-year rate
    that compounds to the realized total return over the sample.
    """
    n = len(returns)
    if n < 1:
        raise ValueError("need at least one return")
    growth = 1.0
    for r in returns:
        growth *= (1.0 + r)
    return growth ** (periods_per_year / n) - 1.0


def annualized_volatility(returns: Sequence[float], periods_per_year=252) -> float:
    """Annualized volatility: the sample standard deviation times
    ``sqrt(periods_per_year)``."""
    if len(returns) < 2:
        raise ValueError("need at least two returns")
    return _std(returns) * math.sqrt(periods_per_year)


def sharpe_ratio(returns: Sequence[float], risk_free=0.0,
                 periods_per_year=252) -> float:
    """Annualized Sharpe ratio of a periodic return series.

    ``(mean_excess / stdev) * sqrt(periods_per_year)`` where ``risk_free`` is the
    per-period risk-free return. Sample standard deviation (ddof=1).
    """
    if len(returns) < 2:
        raise ValueError("need at least two returns")
    excess = [r - risk_free for r in returns]
    sd = _std(excess)
    if sd <= 0.0:
        raise ValueError("zero-variance returns")
    return _mean(excess) / sd * math.sqrt(periods_per_year)


def sortino_ratio(returns: Sequence[float], risk_free=0.0, target=0.0,
                  periods_per_year=252) -> float:
    """Annualized Sortino ratio: excess mean over downside deviation.

    Downside deviation uses only returns below ``target`` (root-mean-square of
    the shortfalls, divided by the full sample count -- the standard
    convention). Raises if there is no downside.
    """
    if len(returns) < 2:
        raise ValueError("need at least two returns")
    excess_mean = _mean([r - risk_free for r in returns])
    downside = [min(r - target, 0.0) for r in returns]
    dd = math.sqrt(sum(d * d for d in downside) / len(returns))
    if dd <= 0.0:
        raise ValueError("no downside deviation (all returns >= target)")
    return excess_mean / dd * math.sqrt(periods_per_year)


def max_drawdown(returns: Sequence[float]) -> float:
    """Maximum peak-to-trough drawdown of the cumulative-return curve.

    Compounds the periodic returns into an equity curve and returns the largest
    fractional drop from a running peak, as a non-negative number (0.2 = a 20%
    drawdown). Empty or all-rising series give 0.
    """
    if not returns:
        return 0.0
    equity = 1.0
    peak = 1.0
    mdd = 0.0
    for r in returns:
        equity *= (1.0 + r)
        if equity > peak:
            peak = equity
        drop = (peak - equity) / peak
        if drop > mdd:
            mdd = drop
    return mdd


def drawdown_curve(returns: Sequence[float]) -> list:
    """Per-period underwater curve: fractional drop from the running peak.

    Compounds the returns into an equity curve and returns, for each period, the
    non-negative drawdown ``(peak - equity)/peak`` at that point (0 at a new
    high). The maximum of this curve is :func:`max_drawdown`.
    """
    out = []
    equity = 1.0
    peak = 1.0
    for r in returns:
        equity *= (1.0 + r)
        if equity > peak:
            peak = equity
        out.append((peak - equity) / peak)
    return out


def longest_drawdown_duration(returns: Sequence[float]) -> int:
    """Longest run of consecutive underwater periods (below a prior peak).

    Counts the maximum number of periods between a peak and the point the equity
    curve first recovers to (or exceeds) it. A series that never falls below its
    running peak returns 0.
    """
    equity = 1.0
    peak = 1.0
    longest = 0
    current = 0
    for r in returns:
        equity *= (1.0 + r)
        if equity >= peak:
            peak = equity
            current = 0
        else:
            current += 1
            if current > longest:
                longest = current
    return longest


def ulcer_index(returns: Sequence[float]) -> float:
    """Ulcer index: RMS of the underwater drawdown curve.

    ``sqrt(mean(drawdown_t^2))`` over the :func:`drawdown_curve` -- a downside risk
    measure that penalizes deep and prolonged drawdowns more than shallow ones,
    unlike volatility which treats up and down moves alike. Zero for a series that
    never draws down.
    """
    dd = drawdown_curve(returns)
    if not dd:
        return 0.0
    return math.sqrt(sum(d * d for d in dd) / len(dd))


def pain_index(returns: Sequence[float]) -> float:
    """Pain index: the average depth of the underwater drawdown curve.

    ``mean(drawdown_t)`` -- the mean fractional distance below the running peak.
    A gentler (L1) cousin of the :func:`ulcer_index` (L2).
    """
    dd = drawdown_curve(returns)
    if not dd:
        return 0.0
    return sum(dd) / len(dd)


def ulcer_performance_index(returns: Sequence[float], risk_free=0.0,
                            periods_per_year=252) -> float:
    """Ulcer performance index (Martin ratio): excess return over the Ulcer index.

    ``(annualized_excess_return) / ulcer_index`` -- a return-per-unit-of-drawdown-
    pain ratio, the drawdown analogue of the Sharpe ratio. Higher is better;
    raises if there is no drawdown (infinite ratio).
    """
    ui = ulcer_index(returns)
    if ui <= 0.0:
        raise ValueError("no drawdown; Ulcer performance index is undefined")
    mean_excess = _mean(returns) - risk_free / periods_per_year
    ann_excess = mean_excess * periods_per_year
    return ann_excess / ui


def pain_ratio(returns: Sequence[float], risk_free=0.0,
               periods_per_year=252) -> float:
    """Pain ratio: annualized excess return over the :func:`pain_index`.

    The L1 analogue of the :func:`ulcer_performance_index`. Higher is better;
    raises when there is no drawdown.
    """
    pi = pain_index(returns)
    if pi <= 0.0:
        raise ValueError("no drawdown; pain ratio is undefined")
    ann_excess = (_mean(returns) - risk_free / periods_per_year) * periods_per_year
    return ann_excess / pi


def rolling_sharpe(returns: Sequence[float], window: int, risk_free=0.0,
                   periods_per_year=252) -> list:
    """Annualized Sharpe ratio over each trailing window of ``window`` periods.

    Returns one Sharpe per window position (``len(returns) - window + 1``
    values), each computed by :func:`sharpe_ratio` on that slice. A
    zero-variance window yields ``float('nan')`` rather than raising, so the
    series stays aligned.
    """
    n = len(returns)
    if window < 2:
        raise ValueError("window must be >= 2")
    if window > n:
        raise ValueError("window longer than the series")
    out = []
    for i in range(n - window + 1):
        chunk = returns[i:i + window]
        try:
            out.append(sharpe_ratio(chunk, risk_free, periods_per_year))
        except ValueError:
            out.append(float("nan"))
    return out


def calmar_ratio(returns: Sequence[float], periods_per_year=252) -> float:
    """Calmar ratio: annualized return divided by the maximum drawdown.

    Annualized return is the geometric ``(prod(1+r))^{periods_per_year/n} - 1``.
    Raises if there is no drawdown (undefined ratio).
    """
    if len(returns) < 2:
        raise ValueError("need at least two returns")
    growth = 1.0
    for r in returns:
        growth *= (1.0 + r)
    ann_return = growth ** (periods_per_year / len(returns)) - 1.0
    mdd = max_drawdown(returns)
    if mdd <= 0.0:
        raise ValueError("no drawdown (Calmar undefined)")
    return ann_return / mdd


def _capture(returns, benchmark_returns, up):
    """Geometric capture ratio over the up (or down) benchmark periods."""
    if len(returns) != len(benchmark_returns):
        raise ValueError("series must be equal length")
    pa = 1.0
    pb = 1.0
    count = 0
    for r, b in zip(returns, benchmark_returns):
        if (b > 0.0) if up else (b < 0.0):
            pa *= (1.0 + r)
            pb *= (1.0 + b)
            count += 1
    if count == 0:
        raise ValueError("no up-market periods" if up else "no down-market periods")
    # Geometric mean per period, then the ratio.
    ga = pa ** (1.0 / count) - 1.0
    gb = pb ** (1.0 / count) - 1.0
    if gb == 0.0:
        raise ValueError("benchmark geometric return is zero over the window")
    return ga / gb


def up_capture(returns, benchmark_returns) -> float:
    """Up-capture ratio: the asset's geometric return in up-benchmark periods
    divided by the benchmark's. Above 1 means the asset outpaces the benchmark
    in rising markets."""
    return _capture(returns, benchmark_returns, up=True)


def down_capture(returns, benchmark_returns) -> float:
    """Down-capture ratio: the asset's geometric return in down-benchmark
    periods over the benchmark's. Below 1 means the asset falls less than the
    benchmark in declining markets (good)."""
    return _capture(returns, benchmark_returns, up=False)


def downside_beta(asset_returns, market_returns) -> float:
    """Beta conditioned on down markets: ``Cov / Var`` over periods where the
    market return is negative.

    Measures how much the asset falls with the market on the downside. Computed
    on the subset of periods with ``market < 0`` using the ordinary sample
    covariance and variance. Raises if there are fewer than two down periods.
    """
    if len(asset_returns) != len(market_returns):
        raise ValueError("series must be equal length")
    a = [ar for ar, m in zip(asset_returns, market_returns) if m < 0.0]
    mk = [m for m in market_returns if m < 0.0]
    n = len(mk)
    if n < 2:
        raise ValueError("need at least two down-market periods")
    ma = sum(a) / n
    mm = sum(mk) / n
    cov = sum((x - ma) * (y - mm) for x, y in zip(a, mk)) / n
    var = sum((y - mm) ** 2 for y in mk) / n
    if var <= 0.0:
        raise ValueError("down-market variance must be positive")
    return cov / var


def tracking_error(returns, benchmark_returns, periods_per_year=252) -> float:
    """Annualized tracking error: stdev of the active (excess) return series.

    ``active_t = r_t - b_t``; the sample standard deviation (ddof=1) scaled by
    ``sqrt(periods_per_year)``. Series must be equal length.
    """
    if len(returns) != len(benchmark_returns):
        raise ValueError("series must be equal length")
    if len(returns) < 2:
        raise ValueError("need at least two observations")
    active = [r - b for r, b in zip(returns, benchmark_returns)]
    return _std(active) * math.sqrt(periods_per_year)


def information_ratio(returns, benchmark_returns, periods_per_year=252) -> float:
    """Information ratio: annualized active return over the tracking error.

    ``mean(active) * periods_per_year / tracking_error`` where the tracking
    error is itself annualized, so this equals
    ``mean(active) / stdev(active) * sqrt(periods_per_year)`` -- the Sharpe of
    the active-return series. Raises if the active returns have no variance.
    """
    if len(returns) != len(benchmark_returns):
        raise ValueError("series must be equal length")
    if len(returns) < 2:
        raise ValueError("need at least two observations")
    active = [r - b for r, b in zip(returns, benchmark_returns)]
    sd = _std(active)
    if sd <= 0.0:
        raise ValueError("zero active-return variance")
    return _mean(active) / sd * math.sqrt(periods_per_year)


def probabilistic_sharpe_ratio(returns, benchmark_sr=0.0):
    """Probabilistic Sharpe ratio (Bailey-López de Prado).

    The probability that the true per-period Sharpe ratio exceeds a ``benchmark_sr``
    (also per period), correcting the estimator's standard error for the sample's
    skewness and (excess) kurtosis and the sample length ``n``:

        PSR = Phi( (SR - SR*) sqrt(n - 1)
                   / sqrt(1 - skew*SR + (kurt-1)/4 * SR^2) ),

    with ``SR`` the per-period Sharpe. Above 0.5 when the observed SR beats the
    benchmark; rises with a longer, less-skewed, thinner-tailed track record.
    """
    from .mathfns import norm_cdf
    n = len(returns)
    if n < 2:
        raise ValueError("need at least two observations")
    mean = _mean(returns)
    sd = _std(returns, ddof=1)
    if sd <= 0.0:
        raise ValueError("zero-variance returns")
    sr = mean / sd
    skew = sample_skewness(returns)
    kurt = sample_kurtosis(returns, excess=False)
    denom = 1.0 - skew * sr + (kurt - 1.0) / 4.0 * sr * sr
    if denom <= 0.0:
        raise ValueError("degenerate PSR denominator")
    z = (sr - benchmark_sr) * math.sqrt(n - 1) / math.sqrt(denom)
    return norm_cdf(z)


def deflated_sharpe_ratio(returns, n_trials, sr_variance=None):
    """Deflated Sharpe ratio (Bailey-López de Prado): PSR against a trials-adjusted benchmark.

    When many strategy variants are tested, the best in-sample Sharpe is inflated by
    selection. The DSR is the :func:`probabilistic_sharpe_ratio` evaluated against a
    benchmark equal to the *expected maximum* of ``n_trials`` independent Sharpe
    estimates with cross-trial variance ``sr_variance``:

        SR* = sqrt(sr_variance) * ((1 - gamma) Phi^{-1}(1 - 1/N)
              + gamma Phi^{-1}(1 - 1/(N e)))

    (``gamma`` the Euler-Mascheroni constant). Lower than the plain PSR for
    ``n_trials > 1``, and falling as more trials are tested. ``sr_variance`` defaults
    to the sampling variance ``1/(n-1)`` of a single per-period Sharpe estimate.
    """
    n = len(returns)
    if n < 2:
        raise ValueError("need at least two observations")
    if n_trials < 1:
        raise ValueError("n_trials must be a positive integer")
    if sr_variance is None:
        sr_variance = 1.0 / (n - 1)
    if sr_variance < 0:
        raise ValueError("sr_variance must be non-negative")
    if n_trials == 1:
        sr_star = 0.0
    else:
        gamma = 0.5772156649015329   # Euler-Mascheroni
        e = math.e
        z1 = _norm_ppf(1.0 - 1.0 / n_trials)
        z2 = _norm_ppf(1.0 - 1.0 / (n_trials * e))
        sr_star = math.sqrt(sr_variance) * ((1.0 - gamma) * z1 + gamma * z2)
    return probabilistic_sharpe_ratio(returns, sr_star)


def minimum_track_record_length(returns, benchmark_sr=0.0, confidence=0.95):
    """Minimum track record length for the Sharpe ratio to beat a benchmark.

    The number of observations at which the :func:`probabilistic_sharpe_ratio`
    would reach ``confidence`` that the true SR exceeds ``benchmark_sr``:

        MinTRL = 1 + (1 - skew*SR + (kurt-1)/4 SR^2) (z_conf / (SR - SR*))^2.

    Requires the observed per-period Sharpe to exceed the benchmark. Longer for a
    smaller edge or a more skewed/fat-tailed series.
    """
    n = len(returns)
    if n < 2:
        raise ValueError("need at least two observations")
    mean = _mean(returns)
    sd = _std(returns, ddof=1)
    if sd <= 0.0:
        raise ValueError("zero-variance returns")
    sr = mean / sd
    if sr <= benchmark_sr:
        raise ValueError("observed Sharpe must exceed the benchmark")
    if not (0.5 < confidence < 1.0):
        raise ValueError("confidence must be in (0.5, 1)")
    skew = sample_skewness(returns)
    kurt = sample_kurtosis(returns, excess=False)
    z = _norm_ppf(confidence)
    factor = 1.0 - skew * sr + (kurt - 1.0) / 4.0 * sr * sr
    return 1.0 + factor * (z / (sr - benchmark_sr)) ** 2


def sample_skewness(returns) -> float:
    """Sample skewness (third standardized moment, population convention).

    ``(1/n) sum (x - mean)^3 / sigma^3`` with the population standard deviation
    (ddof=0). Positive means a longer right tail. Raises on a degenerate
    (zero-variance) series.
    """
    n = len(returns)
    if n < 2:
        raise ValueError("need at least two observations")
    m = _mean(returns)
    var = sum((x - m) ** 2 for x in returns) / n
    if var <= 0.0:
        raise ValueError("zero-variance returns")
    s3 = var ** 1.5
    return sum((x - m) ** 3 for x in returns) / n / s3


def sample_kurtosis(returns, excess=True) -> float:
    """Sample kurtosis (fourth standardized moment, population convention).

    ``(1/n) sum (x - mean)^4 / sigma^4``; with ``excess=True`` subtracts 3 so a
    normal distribution reads 0 (fat tails positive). Raises on zero variance.
    """
    n = len(returns)
    if n < 2:
        raise ValueError("need at least two observations")
    m = _mean(returns)
    var = sum((x - m) ** 2 for x in returns) / n
    if var <= 0.0:
        raise ValueError("zero-variance returns")
    k = sum((x - m) ** 4 for x in returns) / n / (var * var)
    return k - 3.0 if excess else k


def _norm_ppf(p):
    """Inverse standard-normal CDF (Acklam's rational approximation)."""
    if not (0.0 < p < 1.0):
        raise ValueError("p must be in (0, 1)")
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def historical_var(returns, confidence=0.95) -> float:
    """Empirical (historical) Value-at-Risk, as a positive loss.

    The ``(1 - confidence)`` percentile of the return distribution, negated to a
    loss. Uses linear-interpolated order statistics -- no distributional
    assumption. A confidence of 0.95 reports the loss the returns exceed 5% of
    the time.
    """
    if len(returns) < 2:
        raise ValueError("need at least two returns")
    q = _percentile(sorted(returns), (1.0 - confidence) * 100.0)
    return -q


def historical_cvar(returns, confidence=0.95) -> float:
    """Empirical conditional VaR (expected shortfall), as a positive loss.

    The average of the returns at or below the historical-VaR threshold, negated
    -- the mean loss in the worst ``1 - confidence`` of periods. Always at least
    the historical VaR. Falls back to the single worst return when the tail
    holds one observation.
    """
    if len(returns) < 2:
        raise ValueError("need at least two returns")
    s = sorted(returns)
    threshold = _percentile(s, (1.0 - confidence) * 100.0)
    tail = [x for x in s if x <= threshold]
    if not tail:
        tail = [s[0]]
    return -sum(tail) / len(tail)


def cornish_fisher_var(returns, confidence=0.95, horizon=1.0) -> float:
    """Cornish-Fisher (skew/kurtosis-adjusted) Value-at-Risk, as a positive loss.

    Expands the standard-normal quantile ``z`` at ``confidence`` with the sample
    skewness ``S`` and excess kurtosis ``K`` of the returns,

        z_cf = z + (z^2-1) S/6 + (z^3-3z) K/24 - (2z^3-5z) S^2/36,

    evaluated at the lower-tail quantile ``z = Phi^{-1}(1-confidence)`` (a
    negative number), then ``VaR = -(mean*horizon + z_cf*sigma*sqrt(horizon))``
    as a positive loss. For a normal series it reduces to the parametric VaR;
    negative skew and fat tails fatten the left tail and push it above the
    Gaussian VaR.
    """
    n = len(returns)
    if n < 2:
        raise ValueError("need at least two returns")
    mu = _mean(returns)
    sd = _std(returns)
    if sd <= 0.0:
        raise ValueError("zero-variance returns")
    sk = sample_skewness(returns)
    ek = sample_kurtosis(returns, excess=True)
    z = _norm_ppf(1.0 - confidence)   # lower-tail quantile (negative)
    z_cf = (z + (z * z - 1.0) * sk / 6.0
            + (z ** 3 - 3.0 * z) * ek / 24.0
            - (2.0 * z ** 3 - 5.0 * z) * sk * sk / 36.0)
    # Loss VaR at horizon: mean scales linearly, deviation with sqrt(horizon).
    return -(mu * horizon + z_cf * sd * math.sqrt(horizon))


def jarque_bera(returns) -> float:
    """Jarque-Bera test statistic for normality of a return series.

    ``JB = n/6 * (skew^2 + excess_kurt^2/4)``, asymptotically chi-squared with 2
    degrees of freedom under normality. Larger values reject normality (the 5%
    critical value is ~5.99). Uses the population skew/kurtosis.
    """
    n = len(returns)
    if n < 2:
        raise ValueError("need at least two observations")
    sk = sample_skewness(returns)
    ek = sample_kurtosis(returns, excess=True)
    return n / 6.0 * (sk * sk + ek * ek / 4.0)


def _percentile(sorted_vals, p):
    """Linear-interpolated percentile ``p`` in [0, 100] of a sorted list."""
    n = len(sorted_vals)
    if n == 1:
        return sorted_vals[0]
    rank = (p / 100.0) * (n - 1)
    lo = int(math.floor(rank))
    hi = min(lo + 1, n - 1)
    frac = rank - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


def omega_ratio(returns: Sequence[float], threshold=0.0) -> float:
    """Omega ratio: probability-weighted gains over losses about a threshold.

    ``sum(max(r - threshold, 0)) / sum(max(threshold - r, 0))`` -- the ratio of
    upside to downside area relative to ``threshold``. Values above 1 mean more
    gain mass than loss mass. Returns ``inf`` when there is no downside; raises
    if there is neither upside nor downside.
    """
    if not returns:
        raise ValueError("need at least one return")
    up = sum(max(r - threshold, 0.0) for r in returns)
    down = sum(max(threshold - r, 0.0) for r in returns)
    if up == 0.0 and down == 0.0:
        raise ValueError("all returns equal the threshold")
    if down == 0.0:
        return float("inf")
    return up / down


def tail_ratio(returns: Sequence[float], pct=5.0) -> float:
    """Tail ratio: the right tail's magnitude over the left tail's.

    ``|percentile(100 - pct)| / |percentile(pct)|`` -- by default the 95th over
    the 5th percentile (in absolute value). Above 1 means the upside tail is
    fatter than the downside. Raises if the lower tail percentile is zero.
    """
    if len(returns) < 2:
        raise ValueError("need at least two returns")
    if not (0.0 < pct < 50.0):
        raise ValueError("pct must be in (0, 50)")
    s = sorted(returns)
    right = abs(_percentile(s, 100.0 - pct))
    left = abs(_percentile(s, pct))
    if left == 0.0:
        raise ValueError("lower-tail percentile is zero")
    return right / left


def hit_rate(returns: Sequence[float]) -> float:
    """Fraction of periods with a strictly positive return."""
    if not returns:
        raise ValueError("need at least one return")
    return sum(1 for r in returns if r > 0.0) / len(returns)


def profit_factor(returns: Sequence[float]) -> float:
    """Gross profits divided by gross losses (absolute).

    Returns ``inf`` when there are no losing periods. Raises if there are no
    gains and no losses.
    """
    gains = sum(r for r in returns if r > 0.0)
    losses = -sum(r for r in returns if r < 0.0)
    if gains == 0.0 and losses == 0.0:
        raise ValueError("no gains or losses")
    if losses == 0.0:
        return float("inf")
    return gains / losses
