"""Differential evolution: derivative-free global optimization.

Nelder-Mead descends to the *nearest* local minimum; differential evolution (Storn-
Price) searches globally. It evolves a population of candidate vectors, each trial
formed by adding a scaled difference of two others to a third (``mutant = a + F(b - c)``)
and mixing it with the target by crossover, keeping the trial only if it lowers the
objective. This explores multimodal, non-convex landscapes where a local method stalls.
Reproducible via a seeded LCG (no external RNG). Pure standard library.
"""


def _lcg(seed):
    state = seed & 0x7FFFFFFF
    def rand():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return state / 0x80000000
    return rand


def differential_evolution(func, bounds, pop_size=None, F=0.8, cr=0.9,
                           max_iter=1000, tol=1e-12, seed=1234567):
    """Minimize ``func`` over box ``bounds`` by differential evolution.

    ``bounds`` is a list of ``(lo, hi)`` per dimension. ``func`` takes a length-``d``
    list and returns a scalar. ``pop_size`` defaults to ``max(15, 10*d)``; ``F`` is the
    differential weight, ``cr`` the crossover rate. Returns a dict with ``x`` (best
    vector), ``fun`` (its objective), ``n_iter`` and ``converged`` (population spread
    below ``tol``). Deterministic for a fixed ``seed``.
    """
    d = len(bounds)
    if d == 0:
        raise ValueError("need at least one dimension")
    if any(lo > hi for lo, hi in bounds):
        raise ValueError("each bound must have lo <= hi")
    rand = _lcg(seed)
    n = pop_size if pop_size is not None else max(15, 10 * d)
    if n < 4:
        raise ValueError("pop_size must be >= 4")

    def clamp(v):
        return [min(max(v[j], bounds[j][0]), bounds[j][1]) for j in range(d)]

    # Initialize population uniformly in the box.
    pop = [[bounds[j][0] + rand() * (bounds[j][1] - bounds[j][0]) for j in range(d)]
           for _ in range(n)]
    fit = [func(ind) for ind in pop]
    best_i = min(range(n), key=lambda i: fit[i])
    best = list(pop[best_i])
    best_f = fit[best_i]

    n_iter = 0
    for it in range(max_iter):
        n_iter = it + 1
        for i in range(n):
            # Pick three distinct others != i.
            a = b = c = i
            while a == i:
                a = int(rand() * n) % n
            while b == i or b == a:
                b = int(rand() * n) % n
            while c == i or c == a or c == b:
                c = int(rand() * n) % n
            # Mutation + binomial crossover.
            jrand = int(rand() * d) % d
            trial = list(pop[i])
            for j in range(d):
                if rand() < cr or j == jrand:
                    trial[j] = pop[a][j] + F * (pop[b][j] - pop[c][j])
            trial = clamp(trial)
            ft = func(trial)
            if ft <= fit[i]:
                pop[i] = trial
                fit[i] = ft
                if ft < best_f:
                    best_f = ft
                    best = list(trial)
        # Convergence: population objective spread collapsed.
        if max(fit) - min(fit) < tol:
            return {"x": best, "fun": best_f, "n_iter": n_iter, "converged": True}
    return {"x": best, "fun": best_f, "n_iter": n_iter, "converged": False}
