"""2-satisfiability: solving boolean formulas with two literals per clause.

A 2-SAT instance is a conjunction of clauses ``(a OR b)`` over boolean variables. Each
clause ``(a OR b)`` is equivalent to the two implications ``(NOT a -> b)`` and
``(NOT b -> a)``, so the whole formula becomes an implication graph on the ``2n`` literals.
The formula is satisfiable exactly when no variable shares a strongly connected component
with its own negation; when it is, a satisfying assignment reads off the SCC order (a
literal is true when its component is *later* in reverse-topological order than its
negation's). Linear time in the number of variables plus clauses. Pure standard library.
"""

from .scc import strongly_connected_components


class TwoSat:
    """A 2-SAT instance over ``n`` boolean variables ``0 .. n-1``.

    Add clauses with `add_clause` (each literal is a variable index for the positive form or
    its bitwise complement ``~v`` for the negation) or the readable `add_or`. `solve` returns
    a satisfying list of booleans, or ``None`` if the formula is unsatisfiable.
    """

    __slots__ = ("n", "_adj")

    def __init__(self, n):
        if n < 0:
            raise ValueError("n must be non-negative")
        self.n = n
        # 2n nodes: variable v's true-literal is 2v, its false-literal is 2v+1
        self._adj = {i: [] for i in range(2 * n)}

    def _node(self, literal):
        # literal >= 0 -> variable `literal` positive; literal < 0 -> variable ~literal negated
        if literal >= 0:
            v, positive = literal, True
        else:
            v, positive = ~literal, False
        if not 0 <= v < self.n:
            raise IndexError("literal refers to a variable out of range")
        return 2 * v + (0 if positive else 1)

    @staticmethod
    def _neg(node):
        return node ^ 1

    def add_or(self, literal_a, literal_b):
        """Add the clause ``(a OR b)``. Literals are ``v`` (positive) or ``~v`` (negated)."""
        a = self._node(literal_a)
        b = self._node(literal_b)
        # (a OR b) == (NOT a -> b) and (NOT b -> a)
        self._adj[self._neg(a)].append(b)
        self._adj[self._neg(b)].append(a)
        return self

    # alias
    add_clause = add_or

    def add_implication(self, literal_a, literal_b):
        """Add ``a -> b`` (equivalent to the clause ``(NOT a OR b)``)."""
        return self.add_or(~literal_a if literal_a >= 0 else ~literal_a, literal_b)

    def force_true(self, literal):
        """Constrain ``literal`` to be true (the unit clause ``(literal OR literal)``)."""
        return self.add_or(literal, literal)

    def solve(self):
        """Return a satisfying ``[bool] * n`` assignment, or ``None`` if unsatisfiable."""
        comps = strongly_connected_components(self._adj)
        comp_id = {}
        for cid, comp in enumerate(comps):
            for node in comp:
                comp_id[node] = cid
        assignment = [False] * self.n
        for v in range(self.n):
            pos = 2 * v
            neg = 2 * v + 1
            if comp_id[pos] == comp_id[neg]:
                return None                       # variable and its negation collide
            # components come out in reverse topological order, so the smaller id is "later"
            # in topological order; a literal is true when its component precedes its negation
            assignment[v] = comp_id[pos] < comp_id[neg]
        return assignment

    def is_satisfiable(self):
        """True if the formula has a satisfying assignment."""
        return self.solve() is not None
