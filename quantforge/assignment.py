"""Combinatorial optimization: Hungarian assignment and 0/1 knapsack.

Two classic discrete optimizers. The Hungarian algorithm finds the minimum-cost perfect
matching between agents and tasks (assignment problem) in ``O(n^3)``. The 0/1 knapsack
dynamic program maximizes value under a weight budget with each item taken at most once.
Both return the optimum exactly. Pure standard library.
"""


def hungarian(cost):
    """Minimum-cost assignment of ``n`` rows to ``n`` columns (Hungarian algorithm).

    ``cost`` is an ``n x n`` matrix. Returns ``(assignment, total_cost)`` where
    ``assignment[i]`` is the column assigned to row ``i`` and each column is used once.
    Uses the ``O(n^3)`` potentials (Kuhn-Munkres) formulation. For ``maximize``, negate
    the costs before calling.
    """
    n = len(cost)
    if n == 0 or any(len(row) != n for row in cost):
        raise ValueError("cost must be a non-empty square matrix")
    INF = float("inf")
    # 1-indexed potentials method (Jonker-style shortest augmenting path).
    u = [0.0] * (n + 1)
    v = [0.0] * (n + 1)
    p = [0] * (n + 1)          # p[j] = row assigned to column j (0 = unassigned)
    way = [0] * (n + 1)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [INF] * (n + 1)
        used = [False] * (n + 1)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = INF
            j1 = -1
            for j in range(1, n + 1):
                if not used[j]:
                    cur = cost[i0 - 1][j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]
                        j1 = j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while j0:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
    assignment = [0] * n
    for j in range(1, n + 1):
        if p[j]:
            assignment[p[j] - 1] = j - 1
    total = sum(cost[i][assignment[i]] for i in range(n))
    return assignment, total


def knapsack_01(weights, values, capacity):
    """0/1 knapsack: maximize total value with total weight ``<= capacity``.

    Each item (``weights[i]``, ``values[i]``) is taken at most once. Returns
    ``(best_value, chosen_indices)`` via the standard ``O(n * capacity)`` DP. ``capacity``
    and all weights must be non-negative integers.
    """
    n = len(weights)
    if n != len(values):
        raise ValueError("weights and values must have equal length")
    if capacity < 0:
        raise ValueError("capacity must be non-negative")
    # dp[w] = best value at weight budget w; keep a choice trace.
    dp = [0] * (capacity + 1)
    take = [[False] * (capacity + 1) for _ in range(n)]
    for i in range(n):
        wi = weights[i]
        vi = values[i]
        if wi < 0:
            raise ValueError("weights must be non-negative")
        for w in range(capacity, wi - 1, -1):
            if dp[w - wi] + vi > dp[w]:
                dp[w] = dp[w - wi] + vi
                take[i][w] = True
    # Backtrack the chosen items.
    chosen = []
    w = capacity
    for i in range(n - 1, -1, -1):
        if take[i][w]:
            chosen.append(i)
            w -= weights[i]
    chosen.reverse()
    return dp[capacity], chosen
