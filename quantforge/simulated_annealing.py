"""Simulated annealing: global optimization by the Metropolis criterion.

Where differential evolution explores with a population, simulated annealing follows a
single point that sometimes moves *uphill* -- accepting a worse candidate with
probability ``exp(-delta / T)`` -- so it can climb out of local minima. The temperature
``T`` cools geometrically, so early on the walk roams freely and late on it settles into
the basin it has found. A robust, memoryless global search for rugged objectives.
Reproducible via a seeded LCG. Pure standard library.
"""

import math


def _lcg(seed):
    state = seed & 0x7FFFFFFF
    def rand():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return state / 0x80000000
    return rand


def simulated_annealing(func, x0, bounds=None, T0=1.0, cooling=0.995,
                        step=1.0, max_iter=10000, seed=1234567):
    """Minimize ``func`` from start ``x0`` by simulated annealing.

    ``func`` takes a length-``d`` list and returns a scalar. ``bounds`` is an optional
    list of ``(lo, hi)`` per dimension (proposals are clamped into the box). ``T0`` is
    the initial temperature, ``cooling`` the per-iteration geometric decay, ``step`` the
    proposal scale (shrinks with temperature). Returns a dict with ``x`` (best point),
    ``fun`` (its value), ``n_iter`` and ``final_temp``. Deterministic for a fixed seed.
    """
    d = len(x0)
    if d == 0:
        raise ValueError("need at least one dimension")
    if not (0.0 < cooling < 1.0):
        raise ValueError("cooling must be in (0, 1)")
    rand = _lcg(seed)

    def clamp(v):
        if bounds is None:
            return v
        return [min(max(v[j], bounds[j][0]), bounds[j][1]) for j in range(d)]

    cur = clamp([float(x) for x in x0])
    cur_f = func(cur)
    best = list(cur)
    best_f = cur_f
    T = T0

    for it in range(1, max_iter + 1):
        scale = step * (T / T0) ** 0.5
        # Propose a neighbour: perturb each coordinate (Box-Muller-free jitter).
        prop = [cur[j] + scale * (2.0 * rand() - 1.0) for j in range(d)]
        prop = clamp(prop)
        pf = func(prop)
        df = pf - cur_f
        if df <= 0.0 or (T > 0.0 and rand() < math.exp(-df / T)):
            cur, cur_f = prop, pf
            if pf < best_f:
                best, best_f = list(prop), pf
        T *= cooling
    return {"x": best, "fun": best_f, "n_iter": max_iter, "final_temp": T}
