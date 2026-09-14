"""Cycle detection in a functional graph (iterated map ``x -> f(x)``).

Iterating a function from a start value eventually repeats, tracing a "rho" shape: a tail of
length ``mu`` leading into a cycle of length ``lam``. Floyd's tortoise-and-hare and Brent's
algorithm both find ``(mu, lam)`` in ``O(mu + lam)`` time and ``O(1)`` memory -- no set of
visited states -- which is what makes Pollard's rho factoring and PRNG-period analysis
practical. Both are provided, along with a helper returning the actual cycle. Pure standard
library.
"""


def floyd_cycle(f, x0, max_iter=10_000_000):
    """Return ``(mu, lam)`` for the iteration ``x -> f(x)`` from ``x0`` (Floyd's algorithm).

    ``mu`` is the index of the first element on the cycle (tail length); ``lam`` is the cycle
    length. Uses two pointers at speed 1 and 2. Raises if no cycle is found within
    ``max_iter`` steps (only possible for an unbounded state space).
    """
    tortoise = f(x0)
    hare = f(f(x0))
    steps = 0
    while tortoise != hare:
        tortoise = f(tortoise)
        hare = f(f(hare))
        steps += 1
        if steps > max_iter:
            raise ValueError("no cycle found within max_iter")
    # find mu: advance from x0 and meeting point in lockstep
    mu = 0
    tortoise = x0
    while tortoise != hare:
        tortoise = f(tortoise)
        hare = f(hare)
        mu += 1
    # find lam: advance one pointer around the cycle
    lam = 1
    hare = f(tortoise)
    while tortoise != hare:
        hare = f(hare)
        lam += 1
    return mu, lam


def brent_cycle(f, x0, max_iter=10_000_000):
    """Return ``(mu, lam)`` for ``x -> f(x)`` from ``x0`` (Brent's algorithm).

    Same result as :func:`floyd_cycle` but typically fewer function evaluations: it compares
    against a checkpoint whose distance doubles each phase to find ``lam`` first, then ``mu``.
    """
    power = lam = 1
    tortoise = x0
    hare = f(x0)
    steps = 0
    while tortoise != hare:
        if power == lam:
            tortoise = hare
            power *= 2
            lam = 0
        hare = f(hare)
        lam += 1
        steps += 1
        if steps > max_iter:
            raise ValueError("no cycle found within max_iter")
    # find mu: hare leads tortoise by lam, then advance together
    tortoise = hare = x0
    for _ in range(lam):
        hare = f(hare)
    mu = 0
    while tortoise != hare:
        tortoise = f(tortoise)
        hare = f(hare)
        mu += 1
    return mu, lam


def cycle_elements(f, x0):
    """Return the list of states forming the cycle reached by iterating ``f`` from ``x0``."""
    mu, lam = brent_cycle(f, x0)
    x = x0
    for _ in range(mu):
        x = f(x)
    cycle = [x]
    y = f(x)
    while y != x:
        cycle.append(y)
        y = f(y)
    return cycle
