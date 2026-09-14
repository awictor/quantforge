"""Lowest common ancestor on a rooted tree, via binary lifting.

The lowest common ancestor (LCA) of two nodes is the deepest node that is an ancestor of
both. Binary lifting precomputes, for every node, its ancestor ``2^j`` levels up
(``O(n log n)`` build), so an LCA query first lifts the deeper node to its partner's depth
then lifts both in powers of two until they meet -- ``O(log n)`` per query. With depths in
hand, the number of edges between any two nodes follows immediately. Pure standard library.
"""


class LCA:
    """Lowest-common-ancestor index over a rooted tree.

    Construct from an adjacency map ``{node: [neighbors]}`` (an undirected tree) and a root.
    ``query(u, v)`` returns their LCA, ``depth(u)`` the edge count from the root, and
    ``distance(u, v)`` the number of edges on the path between two nodes. Nodes may be any
    hashable label.
    """

    __slots__ = ("_root", "_depth", "_up", "_log", "_index", "_nodes")

    def __init__(self, adjacency, root):
        if root not in adjacency:
            raise ValueError("root is not a node of the tree")
        # assign each label a dense integer index
        self._nodes = list(adjacency)
        self._index = {node: i for i, node in enumerate(self._nodes)}
        n = len(self._nodes)
        self._root = root

        # iterative DFS to get parent and depth (recursion would overflow on deep trees)
        parent = [-1] * n
        depth = [0] * n
        visited = [False] * n
        ri = self._index[root]
        visited[ri] = True
        stack = [ri]
        order = []
        while stack:
            u = stack.pop()
            order.append(u)
            for w in adjacency[self._nodes[u]]:
                wi = self._index[w]
                if not visited[wi]:
                    visited[wi] = True
                    parent[wi] = u
                    depth[wi] = depth[u] + 1
                    stack.append(wi)
        if not all(visited):
            raise ValueError("adjacency is not a single connected tree from root")

        self._depth = depth
        # binary-lifting table: up[j][v] = 2^j-th ancestor of v (-1 past the root)
        self._log = max(1, (n - 1).bit_length())
        up = [[-1] * n for _ in range(self._log)]
        up[0] = parent[:]
        for j in range(1, self._log):
            upj = up[j]
            upjm = up[j - 1]
            for v in range(n):
                mid = upjm[v]
                upj[v] = upjm[mid] if mid != -1 else -1
        self._up = up

    def depth(self, node):
        """Edge count from the root to ``node``."""
        return self._depth[self._index[node]]

    def query(self, u, v):
        """Return the lowest common ancestor of ``u`` and ``v``."""
        iu = self._index[u]
        iv = self._index[v]
        du, dv = self._depth[iu], self._depth[iv]
        # lift the deeper node up to the shallower node's depth
        if du < dv:
            iu, iv = iv, iu
            du, dv = dv, du
        diff = du - dv
        for j in range(self._log):
            if (diff >> j) & 1:
                iu = self._up[j][iu]
        if iu == iv:
            return self._nodes[iu]
        # lift both together until their parents coincide
        for j in range(self._log - 1, -1, -1):
            if self._up[j][iu] != self._up[j][iv]:
                iu = self._up[j][iu]
                iv = self._up[j][iv]
        return self._nodes[self._up[0][iu]]

    def distance(self, u, v):
        """Number of edges on the path between ``u`` and ``v``."""
        w = self.query(u, v)
        return self.depth(u) + self.depth(v) - 2 * self.depth(w)

    def is_ancestor(self, u, v):
        """True if ``u`` is an ancestor of ``v`` (a node is its own ancestor)."""
        return self.query(u, v) == u
