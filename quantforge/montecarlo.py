"""Monte Carlo engine for options, pure standard library.

Simulates geometric Brownian motion with the standard-library Mersenne Twister
(``random.Random``) and offers two variance-reduction techniques:

  * antithetic variates: for each normal draw Z, also use -Z, halving the
    number of independent normals and cancelling odd-moment noise;
  * control variates: subtract a correlated payoff whose expectation is known
    in closed form, then add that expectation back. The arithmetic-average
    Asian uses the geometric-average Asian (Kemna-Vorst) as its control, which
    cuts the standard error by an order of magnitude.

Every price is returned with a Monte Carlo standard error so callers can size
their sample. Passing a ``seed`` makes runs reproducible for testing.
"""

import math
import random
from dataclasses import dataclass

from .bsm import OptionType, _coerce_type, _validate
from .exotics import geometric_asian


@dataclass(frozen=True)
class MCResult:
    price: float
    std_error: float
    n_paths: int

    def confidence_interval(self, z: float = 1.96):
        """Return a (low, high) CI; default z=1.96 is ~95%."""
        half = z * self.std_error
        return (self.price - half, self.price + half)


def _terminal_payoffs(S, K, t, r, sigma, ot, b, n_paths, rng, antithetic):
    """Yield discounted terminal payoffs of a European option (one time step)."""
    drift = (b - 0.5 * sigma * sigma) * t
    vol = sigma * math.sqrt(t)
    disc = math.exp(-r * t)
    sign = 1.0 if ot is OptionType.CALL else -1.0

    def payoff(z):
        sT = S * math.exp(drift + vol * z)
        return disc * max(sign * (sT - K), 0.0)

    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        z = rng.gauss(0.0, 1.0)
        if antithetic:
            yield 0.5 * (payoff(z) + payoff(-z))
        else:
            yield payoff(z)


def _summarize(samples):
    n = len(samples)
    mean = sum(samples) / n
    if n < 2:
        return mean, 0.0
    var = sum((x - mean) ** 2 for x in samples) / (n - 1)
    return mean, math.sqrt(var / n)


def european_mc(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                n_paths=100_000, antithetic=True, seed=None) -> MCResult:
    """Monte Carlo price of a European option (converges to the BSM value)."""
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    rng = random.Random(seed)
    samples = list(_terminal_payoffs(S, K, t, r, sigma, ot, b, n_paths, rng, antithetic))
    price, se = _summarize(samples)
    return MCResult(price=price, std_error=se, n_paths=len(samples))


def _simulate_average_paths(S, t, r, sigma, b, n_steps, n_paths, rng, antithetic):
    """Generate (arithmetic_avg, geometric_avg) of the price path for each run.

    Averages exclude the initial spot and include the value at each of the
    ``n_steps`` monitoring dates up to expiry.
    """
    dt = t / n_steps
    drift = (b - 0.5 * sigma * sigma) * dt
    vol = sigma * math.sqrt(dt)

    def one_path(zs):
        s = S
        arith_sum = 0.0
        log_sum = 0.0
        for z in zs:
            s *= math.exp(drift + vol * z)
            arith_sum += s
            log_sum += math.log(s)
        return arith_sum / n_steps, math.exp(log_sum / n_steps)

    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        zs = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        yield one_path(zs)
        if antithetic:
            yield one_path([-z for z in zs])


def arithmetic_asian_mc(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                        n_steps=50, n_paths=50_000, antithetic=True,
                        control_variate=True, seed=None) -> MCResult:
    """Price a fixed-strike arithmetic-average-price Asian option.

    With ``control_variate=True`` the geometric-average Asian (known in closed
    form) is used as a control, dramatically reducing the standard error since
    the two averages are almost perfectly correlated.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    rng = random.Random(seed)
    disc = math.exp(-r * t)
    sign = 1.0 if ot is OptionType.CALL else -1.0

    geo_closed = geometric_asian(S, K, t, r, sigma, ot, b) if control_variate else 0.0

    arith_payoffs = []
    controlled = []
    for a_avg, g_avg in _simulate_average_paths(
            S, t, r, sigma, b, n_steps, n_paths, rng, antithetic):
        arith_p = disc * max(sign * (a_avg - K), 0.0)
        arith_payoffs.append(arith_p)
        if control_variate:
            geo_p = disc * max(sign * (g_avg - K), 0.0)
            # beta = 1 is near-optimal here (correlation ~ 1); add back the
            # known geometric expectation.
            controlled.append(arith_p - geo_p + geo_closed)

    samples = controlled if control_variate else arith_payoffs
    price, se = _summarize(samples)
    return MCResult(price=price, std_error=se, n_paths=len(samples))
