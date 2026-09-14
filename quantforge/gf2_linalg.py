"""Linear algebra over GF(2): solving XOR/boolean systems by Gaussian elimination.

A linear system modulo 2 -- each equation an XOR of a subset of boolean variables equal to
0 or 1 -- is solved by the same row reduction as over the reals, but with XOR for addition
and no division needed. Rows are bit-packed into Python integers, so elimination is a fast
XOR loop. This underlies Gaussian-elimination attacks, lights-out puzzles, linear codes,
and cycle-space computations. Provides a solver, rank, and a nullspace basis. Pure standard
library.
"""


def gf2_rank(rows):
    """Rank over GF(2) of a matrix given as a list of bit-packed integer rows."""
    basis = []
    for r in rows:
        cur = r
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
            basis.sort(reverse=True)
    return len(basis)


def solve_gf2(equations, rhs, n_vars):
    """Solve an XOR linear system over GF(2). Returns one solution list, or ``None``.

    ``equations[i]`` is a bit-packed integer whose bit ``j`` (value ``1 << j``) marks that
    variable ``j`` appears in equation ``i``; ``rhs[i]`` is that equation's right-hand side
    (0 or 1). ``n_vars`` is the number of variables. Returns a list of ``n_vars`` bits
    (free variables set to 0), or ``None`` if the system is inconsistent.
    """
    if len(equations) != len(rhs):
        raise ValueError("equations and rhs must have equal length")
    # augmented rows: keep the coefficient mask and the rhs bit together
    rows = [(equations[i], rhs[i] & 1) for i in range(len(equations))]
    pivot_row_for_col = {}
    used = [False] * len(rows)
    # forward elimination, one pivot column at a time
    for col in range(n_vars):
        bit = 1 << col
        # find an unused row with this column set
        pr = -1
        for i in range(len(rows)):
            if not used[i] and (rows[i][0] & bit):
                pr = i
                break
        if pr == -1:
            continue
        used[pr] = True
        pivot_row_for_col[col] = pr
        pmask, prhs = rows[pr]
        for i in range(len(rows)):
            if i != pr and (rows[i][0] & bit):
                rows[i] = (rows[i][0] ^ pmask, rows[i][1] ^ prhs)
    # consistency: any row with empty mask but rhs 1 is 0 == 1
    for i in range(len(rows)):
        if rows[i][0] == 0 and rows[i][1] == 1:
            return None
    # back-substitute: pivot columns take their row's rhs, free columns stay 0
    solution = [0] * n_vars
    for col, pr in pivot_row_for_col.items():
        solution[col] = rows[pr][1]
    return solution


def gf2_nullspace_basis(equations, n_vars):
    """Return a basis for the nullspace ``{x : A x = 0}`` over GF(2).

    Each basis vector is a list of ``n_vars`` bits. The number of vectors is
    ``n_vars - rank``. An empty list means only the zero vector satisfies the system.
    """
    rows = [e for e in equations]
    pivot_col_of_row = []
    pivot_cols = set()
    # reduced row echelon form
    reduced = []
    for r in rows:
        cur = r
        for prow, pcol in zip(reduced, pivot_col_of_row):
            if cur & (1 << pcol):
                cur ^= prow
        if cur:
            # pick its lowest set bit as pivot column
            pcol = (cur & -cur).bit_length() - 1
            # eliminate this pivot from earlier rows
            for idx in range(len(reduced)):
                if reduced[idx] & (1 << pcol):
                    reduced[idx] ^= cur
            reduced.append(cur)
            pivot_col_of_row.append(pcol)
            pivot_cols.add(pcol)
    free_cols = [c for c in range(n_vars) if c not in pivot_cols]
    basis = []
    for free in free_cols:
        vec = [0] * n_vars
        vec[free] = 1
        # each pivot row expresses its pivot var in terms of free vars
        for prow, pcol in zip(reduced, pivot_col_of_row):
            if prow & (1 << free):
                vec[pcol] = 1
        basis.append(vec)
    return basis
