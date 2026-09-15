"""Sparse-matrix reordering: reverse Cuthill-McKee bandwidth reduction.

The *bandwidth* of a matrix -- the largest ``|i - j|`` over nonzero ``A[i][j]`` -- controls the
cost and fill-in of banded/sparse factorizations: a small bandwidth means a cheap solve. A
symmetric matrix's bandwidth depends on how its rows/columns are ordered, and the reverse
Cuthill-McKee (RCM) algorithm finds a permutation that dramatically shrinks it. RCM does a
breadth-first traversal of the matrix's adjacency graph starting from a low-degree node, visiting
neighbours in increasing-degree order, then *reverses* the result (George's observation that the
reversal reduces the profile further).

:func:`reverse_cuthill_mckee` returns the permutation; :func:`apply_permutation` symmetrically
permutes the matrix; :func:`matrix_bandwidth` measures the result. Pure standard library.
"""


def matrix_bandwidth(A, tol=0.0):
    """Bandwidth of ``A``: the largest ``|i - j|`` among entries with ``|A[i][j]| > tol``."""
    n = len(A)
    bw = 0
    for i in range(n):
        for j in range(n):
            if abs(A[i][j]) > tol and abs(i - j) > bw:
                bw = abs(i - j)
    return bw


def _adjacency(A, tol=0.0):
    n = len(A)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (abs(A[i][j]) > tol or abs(A[j][i]) > tol):
                adj[i].append(j)
    for i in range(n):
        adj[i] = sorted(set(adj[i]))
    return adj


def reverse_cuthill_mckee(A, tol=0.0):
    """Reverse Cuthill-McKee ordering of the symmetric structure of ``A``.

    Returns a permutation ``perm`` (a list of the original indices in their new order) that,
    applied symmetrically, reduces the bandwidth. Handles disconnected structures by restarting
    from the lowest-degree unvisited node.
    """
    n = len(A)
    adj = _adjacency(A, tol)
    degree = [len(adj[i]) for i in range(n)]
    visited = [False] * n
    order = []
    # process connected components, each seeded by its lowest-degree unvisited node
    remaining = sorted(range(n), key=lambda i: degree[i])
    for seed in remaining:
        if visited[seed]:
            continue
        visited[seed] = True
        queue = [seed]
        head = 0
        while head < len(queue):
            node = queue[head]
            head += 1
            order.append(node)
            # neighbours in increasing degree order
            nbrs = sorted((v for v in adj[node] if not visited[v]), key=lambda v: degree[v])
            for v in nbrs:
                visited[v] = True
                queue.append(v)
    order.reverse()                     # the "reverse" in RCM
    return order


def apply_permutation(A, perm):
    """Symmetrically permute ``A`` by ``perm``: ``B[i][j] = A[perm[i]][perm[j]]``."""
    n = len(perm)
    return [[A[perm[i]][perm[j]] for j in range(n)] for i in range(n)]


def permute_vector(v, perm):
    """Permute a vector: ``w[i] = v[perm[i]]`` (apply the same reordering to a RHS)."""
    return [v[perm[i]] for i in range(len(perm))]


def inverse_permutation(perm):
    """Inverse of a permutation, for mapping a permuted solution back to the original order."""
    inv = [0] * len(perm)
    for new, old in enumerate(perm):
        inv[old] = new
    return inv
