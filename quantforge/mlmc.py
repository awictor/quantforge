"""Multi-level Monte Carlo (Giles 2008) for path-dependent options.

Standard Monte Carlo with a fine time grid pays for accuracy at every path.
Multi-level Monte Carlo instead writes the fine-grid expectation as a telescoping
sum over levels ``l = 0..L`` with step ``dt_l = T / M^l``:

    E[P_L] = E[P_0] + sum_{l=1}^L E[P_l - P_{l-1}],

and estimates each correction with *coupled* fine/coarse paths driven by the
same Brownian increments (the coarse step sums each pair of fine increments).
Because ``P_l - P_{l-1}`` has small variance on fine levels, far fewer samples
are needed there, so the total cost to reach a target accuracy is much lower
than single-level MC. This module applies it to a fixed-strike arithmetic Asian
option and reports the estimate with its per-level sample counts. Pure standard
library.
"""

import math
import random

from .bsm import OptionType, _coerce_type, _validate
from .montecarlo import _summarize, MCResult


def _asian_level(S, K, t, r, sigma, b, sign, level, M, n_paths, rng):
    """Mean and variance of the level-``l`` correction ``P_l - P_{l-1}``.

    At level 0 returns the plain estimator on the coarsest grid; for ``l >= 1``
    the fine grid has ``M^l`` steps and the coarse ``M^{l-1}``, sharing Brownian
    increments (each coarse increment sums ``M`` fine ones).
    """
    disc = math.exp(-r * t)
    n_fine = M ** level
    dt_f = t / n_fine
    sq_f = math.sqrt(dt_f)
    samples = []

    def payoff(avg):
        return disc * max(sign * (avg - K), 0.0)

    if level == 0:
        for _ in range(n_paths):
            s = S
            tot = 0.0
            for _i in range(n_fine):
                s *= math.exp((b - 0.5 * sigma * sigma) * dt_f
                              + sigma * sq_f * rng.gauss(0.0, 1.0))
                tot += s
            samples.append(payoff(tot / n_fine))
        return samples

    n_coarse = M ** (level - 1)
    dt_c = t / n_coarse
    for _ in range(n_paths):
        sf = S
        sc = S
        tot_f = 0.0
        tot_c = 0.0
        c_incr = 0.0          # accumulated fine increments for the coarse step
        for i in range(n_fine):
            z = rng.gauss(0.0, 1.0)
            dW = sq_f * z
            sf *= math.exp((b - 0.5 * sigma * sigma) * dt_f + sigma * dW)
            tot_f += sf
            c_incr += dW
            if (i + 1) % M == 0:
                sc *= math.exp((b - 0.5 * sigma * sigma) * dt_c + sigma * c_incr)
                tot_c += sc
                c_incr = 0.0
        samples.append(payoff(tot_f / n_fine) - payoff(tot_c / n_coarse))
    return samples


def mlmc_asian(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
               levels=4, M=2, n_paths=20_000, seed=None):
    """Arithmetic-average Asian price by multi-level Monte Carlo.

    Runs levels ``0..levels`` with refinement factor ``M`` (fine grid at the top
    has ``M^levels`` steps). ``n_paths`` is the sample count at level 0; deeper
    levels use fewer (``n_paths // M^l``, floored) since their corrections have
    lower variance. Returns an :class:`~quantforge.MCResult` whose ``n_paths`` is
    the total sample count across levels.

    A flat run at ``levels=0`` is plain single-grid Monte Carlo; increasing
    ``levels`` refines the time discretisation while sharing the cost across
    coarser levels.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if levels < 0 or M < 2:
        raise ValueError("need levels >= 0 and M >= 2")
    sign = 1.0 if ot is OptionType.CALL else -1.0
    rng = random.Random(seed)

    total_mean = 0.0
    total_var = 0.0
    total_samples = 0
    for level in range(levels + 1):
        nl = max(n_paths // (M ** level), 500)
        samples = _asian_level(S, K, t, r, sigma, b, sign, level, M, nl, rng)
        m, se = _summarize(samples)
        total_mean += m
        total_var += se * se     # levels are independent -> variances add
        total_samples += len(samples)
    return MCResult(price=total_mean, std_error=math.sqrt(total_var),
                    n_paths=total_samples)
