"""Realized (historical) volatility estimators from a price series.

All estimators return an *annualized* volatility. The range-based estimators
(Parkinson, Garman-Klass, Rogers-Satchell, Yang-Zhang) use intraday
high/low/open/close bars and are far more efficient than close-to-close for the
same sample size; Yang-Zhang additionally handles overnight jumps and drift.

Inputs are plain sequences of floats — no NumPy. ``periods_per_year`` defaults
to 252 (trading days); pass 252*390 for minute bars, etc.

References: Parkinson (1980), Garman-Klass (1980), Rogers-Satchell (1991),
Yang-Zhang (2000).
"""

import math
from dataclasses import dataclass
from typing import Sequence

_LN2 = math.log(2.0)


def _log_returns(closes: Sequence[float]):
    out = []
    for i in range(1, len(closes)):
        if closes[i] <= 0 or closes[i - 1] <= 0:
            raise ValueError("prices must be positive")
        out.append(math.log(closes[i] / closes[i - 1]))
    return out


def close_to_close(closes: Sequence[float], periods_per_year: int = 252,
                   ddof: int = 1) -> float:
    """Classic close-to-close realized volatility (annualized).

    Uses the sample standard deviation of log returns with ``ddof`` degrees of
    freedom removed (1 = unbiased sample variance).
    """
    rets = _log_returns(closes)
    n = len(rets)
    if n <= ddof:
        raise ValueError("not enough returns for the requested ddof")
    mean = sum(rets) / n
    var = sum((r - mean) ** 2 for r in rets) / (n - ddof)
    return math.sqrt(var * periods_per_year)


def ewma_vol(closes: Sequence[float], lam: float = 0.94,
             periods_per_year: int = 252) -> float:
    """RiskMetrics-style exponentially weighted volatility (annualized).

    Variance_t = lam * Variance_{t-1} + (1 - lam) * r_t^2, seeded with the
    first squared return. ``lam=0.94`` is the RiskMetrics daily default.
    """
    if not (0.0 < lam < 1.0):
        raise ValueError("lam must be in (0, 1)")
    rets = _log_returns(closes)
    if not rets:
        raise ValueError("need at least two prices")
    var = rets[0] ** 2
    for r in rets[1:]:
        var = lam * var + (1.0 - lam) * r * r
    return math.sqrt(var * periods_per_year)


def _validate_ohlc(opens, highs, lows, closes):
    n = len(closes)
    if not (len(opens) == len(highs) == len(lows) == n):
        raise ValueError("OHLC series must be the same length")
    if n < 1:
        raise ValueError("need at least one bar")
    for o, h, l, c in zip(opens, highs, lows, closes):
        if min(o, h, l, c) <= 0:
            raise ValueError("prices must be positive")
        if h < max(o, c) or l > min(o, c):
            raise ValueError("bar violates high/low bounds")
    return n


def parkinson(highs: Sequence[float], lows: Sequence[float],
              periods_per_year: int = 252) -> float:
    """Parkinson high-low range estimator (annualized).

    var = (1 / (4 ln2)) * mean( ln(H/L)^2 ). ~5x more efficient than
    close-to-close but ignores drift and overnight moves.
    """
    n = len(highs)
    if n < 1 or len(lows) != n:
        raise ValueError("need matching high/low series")
    acc = 0.0
    for h, l in zip(highs, lows):
        if h <= 0 or l <= 0 or h < l:
            raise ValueError("invalid high/low bar")
        acc += math.log(h / l) ** 2
    var = acc / (4.0 * _LN2 * n)
    return math.sqrt(var * periods_per_year)


def garman_klass(opens, highs, lows, closes, periods_per_year: int = 252) -> float:
    """Garman-Klass OHLC estimator (annualized).

    var = mean( 0.5*ln(H/L)^2 - (2ln2 - 1)*ln(C/O)^2 ). Uses the full bar; more
    efficient than Parkinson, still assumes no overnight jump or drift.
    """
    n = _validate_ohlc(opens, highs, lows, closes)
    acc = 0.0
    for o, h, l, c in zip(opens, highs, lows, closes):
        hl = math.log(h / l)
        co = math.log(c / o)
        acc += 0.5 * hl * hl - (2.0 * _LN2 - 1.0) * co * co
    var = acc / n
    return math.sqrt(var * periods_per_year)


def rogers_satchell(opens, highs, lows, closes, periods_per_year: int = 252) -> float:
    """Rogers-Satchell OHLC estimator (annualized).

    var = mean( ln(H/C)ln(H/O) + ln(L/C)ln(L/O) ). Drift-independent: stays
    unbiased even when the underlying has a non-zero mean return.
    """
    n = _validate_ohlc(opens, highs, lows, closes)
    acc = 0.0
    for o, h, l, c in zip(opens, highs, lows, closes):
        acc += (math.log(h / c) * math.log(h / o)
                + math.log(l / c) * math.log(l / o))
    var = acc / n
    return math.sqrt(var * periods_per_year)


def yang_zhang(opens, highs, lows, closes, periods_per_year: int = 252) -> float:
    """Yang-Zhang estimator (annualized): drift-independent and jump-robust.

    Combines overnight (close-to-open) variance, open-to-close variance, and
    the Rogers-Satchell term:

        var = var_overnight + k * var_open_to_close + (1 - k) * var_RS
        k   = 0.34 / (1.34 + (N+1)/(N-1))

    Requires the previous close, so bars are chained: overnight return uses
    ln(O_t / C_{t-1}).
    """
    n = _validate_ohlc(opens, highs, lows, closes)
    if n < 2:
        raise ValueError("Yang-Zhang needs at least two bars")

    # Overnight (close-to-open) log returns, chained across bars.
    overnight = [math.log(opens[i] / closes[i - 1]) for i in range(1, n)]
    # Open-to-close log returns on the same bars.
    open_close = [math.log(closes[i] / opens[i]) for i in range(1, n)]
    m = len(overnight)  # = n - 1

    mean_on = sum(overnight) / m
    var_on = sum((x - mean_on) ** 2 for x in overnight) / (m - 1) if m > 1 else 0.0
    mean_oc = sum(open_close) / m
    var_oc = sum((x - mean_oc) ** 2 for x in open_close) / (m - 1) if m > 1 else 0.0

    # Rogers-Satchell over the same bars (indices 1..n-1).
    rs = 0.0
    for i in range(1, n):
        o, h, l, c = opens[i], highs[i], lows[i], closes[i]
        rs += (math.log(h / c) * math.log(h / o)
               + math.log(l / c) * math.log(l / o))
    var_rs = rs / m

    k = 0.34 / (1.34 + (m + 1) / (m - 1)) if m > 1 else 0.34
    var = var_on + k * var_oc + (1.0 - k) * var_rs
    return math.sqrt(var * periods_per_year)


@dataclass(frozen=True)
class VolReport:
    close_to_close: float
    parkinson: float
    garman_klass: float
    rogers_satchell: float
    yang_zhang: float
    ewma: float


def vol_report(opens, highs, lows, closes, periods_per_year: int = 252,
               ewma_lambda: float = 0.94) -> VolReport:
    """Compute every estimator at once for an OHLC series."""
    return VolReport(
        close_to_close=close_to_close(closes, periods_per_year),
        parkinson=parkinson(highs, lows, periods_per_year),
        garman_klass=garman_klass(opens, highs, lows, closes, periods_per_year),
        rogers_satchell=rogers_satchell(opens, highs, lows, closes, periods_per_year),
        yang_zhang=yang_zhang(opens, highs, lows, closes, periods_per_year),
        ewma=ewma_vol(closes, ewma_lambda, periods_per_year),
    )


@dataclass(frozen=True)
class VolConePoint:
    window: int              # rolling window length (in bars)
    minimum: float
    p25: float
    median: float
    p75: float
    maximum: float
    current: float           # most-recent window's realized vol


def _percentile(sorted_vals, p):
    """Linear-interpolated percentile of an already-sorted list (p in [0,1])."""
    n = len(sorted_vals)
    if n == 1:
        return sorted_vals[0]
    idx = p * (n - 1)
    lo = int(math.floor(idx))
    hi = min(lo + 1, n - 1)
    frac = idx - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


def vol_cone(closes: Sequence[float], windows: Sequence[int],
             periods_per_year: int = 252):
    """Realized-volatility cone: the distribution of rolling realized vol per window.

    For each window length, computes the annualized close-to-close realized vol
    over every rolling block of returns of that length, then reports the min,
    25th/50th/75th percentiles, max, and the most-recent (current) value. This
    is the standard "vol cone" used to judge whether current realized vol is
    high or low versus its own history at each horizon.

    Args:
        closes: the price series.
        windows: rolling window lengths in *returns* (e.g. [5, 21, 63, 126]).
        periods_per_year: annualization factor.

    Returns a list of :class:`VolConePoint`, one per window (skipping windows
    too long for the data).
    """
    rets = _log_returns(closes)
    n = len(rets)
    out = []
    for w in sorted(windows):
        if w < 2 or w > n:
            continue
        vols = []
        for start in range(0, n - w + 1):
            block = rets[start:start + w]
            mean = sum(block) / w
            var = sum((x - mean) ** 2 for x in block) / (w - 1)
            vols.append(math.sqrt(var * periods_per_year))
        vols_sorted = sorted(vols)
        out.append(VolConePoint(
            window=w,
            minimum=vols_sorted[0],
            p25=_percentile(vols_sorted, 0.25),
            median=_percentile(vols_sorted, 0.5),
            p75=_percentile(vols_sorted, 0.75),
            maximum=vols_sorted[-1],
            current=vols[-1],   # most-recent rolling block
        ))
    return out


@dataclass(frozen=True)
class GarchParams:
    omega: float             # long-run variance intercept
    alpha: float             # weight on last squared return (ARCH)
    beta: float              # weight on last variance (GARCH)

    @property
    def persistence(self):
        return self.alpha + self.beta

    @property
    def long_run_variance(self):
        p = self.persistence
        return self.omega / (1.0 - p) if p < 1.0 else float("inf")


def fit_garch(returns: Sequence[float], periods_per_year: int = 252):
    """Fit a GARCH(1,1) variance model to a return series by quasi-MLE.

    Model: ``h_t = omega + alpha * r_{t-1}^2 + beta * h_{t-1}`` with the returns
    assumed conditionally normal (Gaussian quasi-likelihood). Fitted with the
    built-in Nelder-Mead over a smooth reparametrization that keeps
    ``omega > 0``, ``alpha, beta >= 0`` and ``alpha + beta < 1`` (stationary).

    Returns ``GarchParams`` (per-period variance parameters).
    """
    from .optimize import nelder_mead

    r = [float(x) for x in returns]
    n = len(r)
    if n < 20:
        raise ValueError("need at least ~20 returns to fit GARCH")
    sample_var = sum((x - sum(r) / n) ** 2 for x in r) / n

    def softplus(x):
        return math.log1p(math.exp(-abs(x))) + max(x, 0.0)

    def unpack(p):
        # omega = softplus; alpha, beta via a softmax-like split of persistence.
        omega = softplus(p[0]) + 1e-12
        # persistence in (0,1) via logistic; split between alpha/beta by logistic.
        pers = 1.0 / (1.0 + math.exp(-p[1]))
        frac = 1.0 / (1.0 + math.exp(-p[2]))
        alpha = pers * frac
        beta = pers * (1.0 - frac)
        return omega, alpha, beta

    def neg_loglik(p):
        omega, alpha, beta = unpack(p)
        h = sample_var
        ll = 0.0
        for i in range(n):
            if h <= 0:
                return 1e18
            ll += 0.5 * (math.log(h) + r[i] * r[i] / h)
            h = omega + alpha * r[i] * r[i] + beta * h
        return ll

    x0 = [math.log(math.expm1(max(sample_var * 0.1, 1e-10))), 0.0, 0.0]
    best, _ = nelder_mead(neg_loglik, x0, step=0.5, max_iter=4000, tol=1e-12)
    omega, alpha, beta = unpack(best)
    return GarchParams(omega=omega, alpha=alpha, beta=beta)


def garch_forecast(params: GarchParams, last_return, last_variance,
                   horizon=1, periods_per_year: int = 252):
    """Forecast annualized volatility ``horizon`` periods ahead under GARCH(1,1).

    The one-step-ahead variance is ``h_1 = omega + alpha r^2 + beta h``. Beyond
    that the expected variance mean-reverts toward the long-run level at rate
    ``persistence`` per step: ``E[h_k] = LR + persistence^{k-1} (h_1 - LR)``.
    Returns the annualized volatility for the ``horizon``-step-ahead period.
    """
    if horizon < 1:
        raise ValueError("horizon must be >= 1")
    h1 = params.omega + params.alpha * last_return * last_return + params.beta * last_variance
    lr = params.long_run_variance
    p = params.persistence
    if horizon == 1 or p >= 1.0:
        h = h1
    else:
        h = lr + (p ** (horizon - 1)) * (h1 - lr)
    return math.sqrt(h * periods_per_year)


def garch_term_variance(params: GarchParams, last_return, last_variance,
                        horizon, periods_per_year: int = 252):
    """Annualized GARCH term (average) volatility over the next ``horizon`` steps.

    An option maturing in ``horizon`` periods is priced off the *average* of the
    per-step conditional variances, not a single step. Summing the mean-reverting
    forecasts ``E[h_k] = LR + persistence^{k-1} (h_1 - LR)`` gives, for
    persistence ``p < 1``,

        avg_var = LR + (h_1 - LR)/horizon * (1 - p^horizon)/(1 - p),

    the closed form of the geometric-series average. This is the volatility to
    feed a Black-Scholes price for that maturity. Returns the annualized term
    volatility ``sqrt(avg_var * periods_per_year)``.
    """
    if horizon < 1:
        raise ValueError("horizon must be >= 1")
    h1 = params.omega + params.alpha * last_return * last_return + params.beta * last_variance
    lr = params.long_run_variance
    p = params.persistence
    if p >= 1.0:
        avg = h1  # unit-root: variance neither reverts nor is finite in LR
    else:
        # sum_{k=1}^{n} E[h_k] = n*LR + (h1-LR) * (1 - p^n)/(1 - p); average it.
        geom = (1.0 - p ** horizon) / (1.0 - p)
        avg = lr + (h1 - lr) * geom / horizon
    return math.sqrt(avg * periods_per_year)


def garch_option_price(params: GarchParams, last_return, last_variance,
                       S, K, r, option_type="call", horizon=None, t=None,
                       periods_per_year: int = 252, b=None):
    """Black-Scholes price using the GARCH term volatility for the maturity.

    Bridges the GARCH variance forecast to an option price: the annualized term
    (average) volatility over ``horizon`` steps -- :func:`garch_term_variance` --
    is fed into the Black-Scholes formula. ``horizon`` is the number of GARCH
    steps to expiry; the option's year fraction ``t`` defaults to
    ``horizon / periods_per_year`` but may be passed explicitly (e.g. to use
    calendar rather than trading time). This lets a fitted GARCH model price
    options consistently with its own vol term structure -- capturing the vol
    mean-reversion that a single spot vol misses.
    """
    from .bsm import price as bsm_price
    if horizon is None:
        raise ValueError("horizon (steps to expiry) is required")
    sigma = garch_term_variance(params, last_return, last_variance, horizon,
                                periods_per_year)
    tt = (horizon / periods_per_year) if t is None else t
    return bsm_price(S, K, tt, r, sigma, option_type, b=b)
