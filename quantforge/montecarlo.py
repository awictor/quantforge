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


def capped_cliquet_mc(S, t, r, sigma, reset_times, local_cap=None,
                      local_floor=0.0, global_cap=None, global_floor=0.0,
                      b=None, n_paths=50_000, antithetic=True, seed=None) -> MCResult:
    """Monte Carlo a locally- and globally-capped cliquet (ratchet).

    The payoff sums the periodic returns of the underlying over consecutive
    reset windows, clipping each period return to ``[local_floor, local_cap]``,
    then clips the running sum to ``[global_floor, global_cap]``. The result is
    discounted at ``r``. This is the standard capped-cliquet structured note;
    the caps make it path-dependent with no closed form.

    Args:
        reset_times: increasing schedule, e.g. [0.25, 0.5, 0.75, 1.0]; the first
            period runs from now (t=0) to reset_times[0].
        local_cap / local_floor: per-period return bounds (cap None = uncapped).
        global_cap / global_floor: bounds on the summed payoff.
    """
    _validate(S, S, t, sigma)
    if b is None:
        b = r
    times = [0.0] + list(reset_times)
    if any(times[i] >= times[i + 1] for i in range(len(times) - 1)):
        raise ValueError("reset_times must be strictly increasing and positive")
    disc = math.exp(-r * t)
    rng = random.Random(seed)

    def clip(x, lo, hi):
        if lo is not None:
            x = max(x, lo)
        if hi is not None:
            x = min(x, hi)
        return x

    def one_path(draws):
        total = 0.0
        for i, z in enumerate(draws):
            dt = times[i + 1] - times[i]
            drift = (b - 0.5 * sigma * sigma) * dt
            ret = math.exp(drift + sigma * math.sqrt(dt) * z) - 1.0
            total += clip(ret, local_floor, local_cap)
        return disc * clip(total, global_floor, global_cap)

    n_periods = len(reset_times)
    samples = []
    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        zs = [rng.gauss(0.0, 1.0) for _ in range(n_periods)]
        samples.append(one_path(zs))
        if antithetic:
            samples.append(one_path([-z for z in zs]))

    price, se = _summarize(samples)
    return MCResult(price=price, std_error=se, n_paths=len(samples))


def barrier_digital_mc(S, K, H, t, r, sigma, option_type=OptionType.CALL,
                       barrier="up-in", b=None, cash=1.0, n_steps=100,
                       n_paths=50_000, antithetic=True, seed=None) -> MCResult:
    """Monte Carlo a cash-or-nothing digital contingent on a barrier condition.

    Pays ``cash`` at expiry if the option finishes in the money (call: S_T > K;
    put: S_T < K) AND the barrier condition holds over the monitored path:

      * ``"up-in"``   / ``"down-in"``   : the barrier H must be touched;
      * ``"up-out"``  / ``"down-out"``  : the barrier H must NOT be touched.

    "up" barriers watch for S >= H, "down" for S <= H. This is the standard
    barrier-contingent binary; the path dependence has no simple closed form.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    barrier = str(barrier).lower()
    if barrier not in ("up-in", "down-in", "up-out", "down-out"):
        raise ValueError("barrier must be up-in/down-in/up-out/down-out")
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")

    up = barrier.startswith("up")
    knock_in = barrier.endswith("in")
    itm = (lambda s: s > K) if ot is OptionType.CALL else (lambda s: s < K)

    dt = t / n_steps
    drift = (b - 0.5 * sigma * sigma) * dt
    vol = sigma * math.sqrt(dt)
    disc = math.exp(-r * t)
    rng = random.Random(seed)

    def one_path(zs):
        s = S
        touched = (up and S >= H) or (not up and S <= H)
        for z in zs:
            s *= math.exp(drift + vol * z)
            if up and s >= H:
                touched = True
            elif not up and s <= H:
                touched = True
        barrier_ok = touched if knock_in else not touched
        return disc * cash if (barrier_ok and itm(s)) else 0.0

    samples = []
    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        zs = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        samples.append(one_path(zs))
        if antithetic:
            samples.append(one_path([-z for z in zs]))

    price, se = _summarize(samples)
    return MCResult(price=price, std_error=se, n_paths=len(samples))


def parisian_barrier_mc(S, K, H, t, r, sigma, window, option_type=OptionType.CALL,
                        barrier="down-out", b=None, n_steps=252, n_paths=40_000,
                        antithetic=True, seed=None) -> MCResult:
    """Monte Carlo a Parisian barrier option.

    Unlike a standard barrier (triggered by a single touch), a Parisian barrier
    triggers only if the spot stays on the barrier's far side for a *consecutive*
    elapsed time of at least ``window`` years. This makes the option robust to
    brief spikes through the level.

    ``barrier`` is one of ``down-out``/``down-in``/``up-out``/``up-in``. "down"
    watches for S <= H, "up" for S >= H. Knock-out pays the vanilla payoff
    unless the barrier is activated; knock-in pays only if it is.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    barrier = str(barrier).lower()
    if barrier not in ("down-out", "down-in", "up-out", "up-in"):
        raise ValueError("barrier must be down-out/down-in/up-out/up-in")
    if window <= 0 or window > t:
        raise ValueError("window must be in (0, t]")
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")

    up = barrier.startswith("up")
    knock_in = barrier.endswith("in")
    dt = t / n_steps
    window_steps = max(1, int(round(window / dt)))
    drift = (b - 0.5 * sigma * sigma) * dt
    vol = sigma * math.sqrt(dt)
    disc = math.exp(-r * t)
    sign = 1.0 if ot is OptionType.CALL else -1.0

    def one_path(zs):
        s = S
        consec = 0
        activated = False
        for z in zs:
            s *= math.exp(drift + vol * z)
            beyond = (s >= H) if up else (s <= H)
            if beyond:
                consec += 1
                if consec >= window_steps:
                    activated = True
            else:
                consec = 0
        payoff = max(sign * (s - K), 0.0)
        alive = activated if knock_in else (not activated)
        return disc * payoff if alive else 0.0

    rng = random.Random(seed)
    samples = []
    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        zs = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        samples.append(one_path(zs))
        if antithetic:
            samples.append(one_path([-z for z in zs]))

    price, se = _summarize(samples)
    return MCResult(price=price, std_error=se, n_paths=len(samples))


def local_vol_mc(S, K, t, r, local_vol_fn, option_type=OptionType.CALL, q=0.0,
                 n_steps=100, n_paths=50_000, antithetic=True, seed=None) -> MCResult:
    """Monte Carlo a European option under a Dupire local-volatility surface.

    Args:
        local_vol_fn: callable ``sigma_loc(S, t_now)`` giving the instantaneous
            local volatility at spot ``S`` and elapsed time ``t_now``.
        q: continuous dividend yield (drift is ``r - q``).

    Evolves ``dS = (r - q) S dt + sigma_loc(S, t) S dW`` with an Euler step in
    log-space. For a flat local vol this reproduces the Black-Scholes price; for
    a genuine Dupire surface the discretely-simulated price is consistent with
    that surface's vanilla prices.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, 0.1)   # sigma checked inside local_vol_fn per step
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    dt = t / n_steps
    sqdt = math.sqrt(dt)
    disc = math.exp(-r * t)
    sign = 1.0 if ot is OptionType.CALL else -1.0
    mu = r - q

    def one_path(zs):
        s = S
        tau = 0.0
        for z in zs:
            sig = local_vol_fn(s, tau)
            if sig < 0:
                raise ValueError("local vol must be non-negative")
            s *= math.exp((mu - 0.5 * sig * sig) * dt + sig * sqdt * z)
            tau += dt
        return disc * max(sign * (s - K), 0.0)

    rng = random.Random(seed)
    samples = []
    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        zs = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        samples.append(one_path(zs))
        if antithetic:
            samples.append(one_path([-z for z in zs]))

    price, se = _summarize(samples)
    return MCResult(price=price, std_error=se, n_paths=len(samples))


def average_strike_asian_mc(S, t, r, sigma, option_type=OptionType.CALL, b=None,
                            n_steps=50, n_paths=50_000, antithetic=True,
                            seed=None) -> MCResult:
    """Monte Carlo an average-strike Asian option.

    The strike is the realized arithmetic average of the monitored path, so a
    call pays ``max(S_T - A, 0)`` and a put ``max(A - S_T, 0)``, where ``A`` is
    the average over the ``n_steps`` monitoring dates. There is no simple closed
    form; the average and terminal spot come from the same simulated path.
    """
    ot = _coerce_type(option_type)
    _validate(S, S, t, sigma)
    if b is None:
        b = r
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    dt = t / n_steps
    drift = (b - 0.5 * sigma * sigma) * dt
    vol = sigma * math.sqrt(dt)
    disc = math.exp(-r * t)
    sign = 1.0 if ot is OptionType.CALL else -1.0

    def one_path(zs):
        s = S
        avg_sum = 0.0
        for z in zs:
            s *= math.exp(drift + vol * z)
            avg_sum += s
        avg = avg_sum / n_steps
        return disc * max(sign * (s - avg), 0.0)   # strike = realized average

    rng = random.Random(seed)
    samples = []
    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        zs = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        samples.append(one_path(zs))
        if antithetic:
            samples.append(one_path([-z for z in zs]))

    price, se = _summarize(samples)
    return MCResult(price=price, std_error=se, n_paths=len(samples))


def autocallable_mc(S, t, r, sigma, observation_times, autocall_barrier,
                    coupon, protection_barrier=None, notional=1.0, b=None,
                    n_paths=50_000, antithetic=True, seed=None) -> MCResult:
    """Monte Carlo an autocallable structured note.

    At each observation date, if the spot is at or above ``autocall_barrier``
    the note redeems early paying ``notional * (1 + coupon * k)`` where ``k`` is
    the observation number (accrued coupons), discounted to today. If it never
    autocalls, at maturity the holder gets the notional back unless the spot
    finished below ``protection_barrier`` (a down-and-in put on the notional),
    in which case they take the downside ``notional * S_T / S``.

    Args:
        observation_times: increasing dates (years); the last is maturity.
        autocall_barrier / protection_barrier: spot levels (same units as S).
        coupon: coupon rate paid per elapsed observation on early redemption.
    """
    _validate(S, S, t, sigma)
    if b is None:
        b = r
    obs = list(observation_times)
    if any(obs[i] >= obs[i + 1] for i in range(len(obs) - 1)) or obs[0] <= 0:
        raise ValueError("observation_times must be strictly increasing and positive")
    if abs(obs[-1] - t) > 1e-9:
        raise ValueError("last observation must be the maturity t")

    steps = [obs[0]] + [obs[i + 1] - obs[i] for i in range(len(obs) - 1)]
    rng = random.Random(seed)

    def one_path(zs):
        s = S
        tau = 0.0
        for k, (dt, z) in enumerate(zip(steps, zs), start=1):
            drift = (b - 0.5 * sigma * sigma) * dt
            s *= math.exp(drift + sigma * math.sqrt(dt) * z)
            tau += dt
            if s >= autocall_barrier and k < len(steps):
                # Early redemption: notional + accrued coupons.
                return math.exp(-r * tau) * notional * (1.0 + coupon * k)
        # Reached maturity without autocalling.
        disc = math.exp(-r * t)
        if protection_barrier is not None and s < protection_barrier:
            return disc * notional * (s / S)          # downside participation
        return disc * notional * (1.0 + coupon * len(steps))

    samples = []
    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        zs = [rng.gauss(0.0, 1.0) for _ in range(len(steps))]
        samples.append(one_path(zs))
        if antithetic:
            samples.append(one_path([-z for z in zs]))

    price, se = _summarize(samples)
    return MCResult(price=price, std_error=se, n_paths=len(samples))


def double_knockout_mc(S, K, t, r, sigma, lower, upper, option_type=OptionType.CALL,
                       b=None, rebate=0.0, n_steps=100, n_paths=50_000,
                       antithetic=True, seed=None) -> MCResult:
    """Monte Carlo a double-knockout barrier option (a corridor).

    The option pays the vanilla payoff only if the spot stays strictly inside
    ``(lower, upper)`` for the whole monitored path; if either barrier is
    breached it knocks out and pays the cash ``rebate`` (at expiry, discounted).
    Also known as a double-barrier knock-out or "corridor" option.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if not (lower < S < upper):
        raise ValueError("require lower < S < upper")
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    dt = t / n_steps
    drift = (b - 0.5 * sigma * sigma) * dt
    vol = sigma * math.sqrt(dt)
    disc = math.exp(-r * t)
    sign = 1.0 if ot is OptionType.CALL else -1.0

    def one_path(zs):
        s = S
        knocked = False
        for z in zs:
            s *= math.exp(drift + vol * z)
            if s <= lower or s >= upper:
                knocked = True
                break
        if knocked:
            return disc * rebate
        return disc * max(sign * (s - K), 0.0)

    rng = random.Random(seed)
    samples = []
    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        zs = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        samples.append(one_path(zs))
        if antithetic:
            samples.append(one_path([-z for z in zs]))

    price, se = _summarize(samples)
    return MCResult(price=price, std_error=se, n_paths=len(samples))
