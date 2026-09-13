"""Series acceleration: Wynn's epsilon algorithm and the Euler transform.

Aitken's delta-squared is a single extrapolation; these push further:

  * ``wynn_epsilon`` -- Wynn's epsilon algorithm, which iterates the Shanks
    transformation to *all* orders in one triangular table. The even-column entries are
    successive Shanks transforms; the last usable one is typically accurate to many
    more digits than Aitken, and it is the standard engine for accelerating slowly- or
    alternating-converging series.
  * ``euler_transform`` -- Euler's transform of an alternating series
    ``sum (-1)^k a_k``, reweighting the terms by forward differences so an alternating
    tail collapses geometrically.

Both take the sequence of *partial sums* (Wynn) or the term sequence (Euler). Pure
standard library.
"""


def wynn_epsilon(partial_sums):
    """Accelerate a sequence of partial sums with Wynn's epsilon algorithm.

    ``partial_sums`` is the list ``[s_0, s_1, ...]``. Returns the best (last stable)
    even-column estimate of the limit. The epsilon table is built with the recurrence
    ``eps[k+1][j] = eps[k-1][j+1] + 1 / (eps[k][j+1] - eps[k][j])``; the even columns
    hold the accelerated limits.
    """
    n = len(partial_sums)
    if n < 3:
        raise ValueError("need at least 3 partial sums")

    # e[k] is the current column; keep previous column for the recurrence.
    prev = [0.0] * (n + 1)         # epsilon_{-1} = 0
    cur = list(partial_sums)       # epsilon_0 = partial sums
    # Collect the last entry of each even column as a limit estimate; once the table
    # degrades (round-off after convergence) later estimates get worse, so return the
    # one whose successive change is smallest rather than blindly the last.
    estimates = [cur[-1]]
    for k in range(1, n):
        nxt = []
        for j in range(n - k):
            denom = cur[j + 1] - cur[j]
            if denom == 0.0:
                nxt.append(prev[j + 1])       # avoid division by zero; carry across
            else:
                nxt.append(prev[j + 1] + 1.0 / denom)
        if k % 2 == 0 and nxt:
            estimates.append(nxt[-1])
        prev = cur
        cur = nxt
        if len(cur) < 2:
            break

    if len(estimates) == 1:
        return estimates[0]
    # Pick the estimate at the step of smallest successive change (best convergence).
    best_idx = min(range(1, len(estimates)),
                   key=lambda i: abs(estimates[i] - estimates[i - 1]))
    return estimates[best_idx]


def euler_transform(terms):
    """Euler transform of an alternating series ``sum_{k>=0} (-1)^k terms[k]``.

    ``terms`` are the non-negative magnitudes ``a_k`` (the alternating sign is applied
    internally). Returns the accelerated estimate of the sum via
    ``sum_k (-1)^k a_k = sum_n (-1)^0 * Delta^n a_0 / 2^{n+1}`` -- forward differences
    reweighted by powers of one half, which converges geometrically even when the raw
    alternating series crawls.
    """
    a = list(terms)
    n = len(a)
    if n == 0:
        raise ValueError("need at least one term")
    # Euler: sum (-1)^k a_k = sum_n (-1)^n Delta^n a_0 / 2^{n+1}, with Delta the
    # forward difference of the magnitudes. Build the difference columns, row 0 gives
    # Delta^n a_0.
    diff = list(a)
    total = 0.0
    half = 0.5
    sign = 1.0
    for _ in range(n):
        total += sign * half * diff[0]
        half *= 0.5
        sign = -sign
        diff = [diff[i + 1] - diff[i] for i in range(len(diff) - 1)]
        if not diff:
            break
    return total
