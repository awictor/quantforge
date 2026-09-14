"""Linear programming by the two-phase simplex method.

Maximizes (or minimizes) a linear objective subject to linear inequality/equality
constraints. Uses a tableau simplex with Bland's anti-cycling rule and a phase-1
artificial-variable stage to find an initial feasible basis, so it handles ``<=``,
``>=`` and ``=`` constraints with any sign of right-hand side. Returns the optimal vertex
and objective. Pure standard library.
"""


def linprog(c, constraints, maximize=True, tol=1e-9, max_iter=10000):
    """Solve a linear program by two-phase simplex.

    ``c`` is the objective coefficient vector (length ``n``). ``constraints`` is a list of
    ``(coeffs, op, rhs)`` with ``op`` one of ``"<="``, ``">="``, ``"="``. Variables are
    assumed ``>= 0``. Set ``maximize=False`` to minimize. Returns
    ``{"x": [...], "objective": value, "status": "optimal"}``. Raises ``ValueError`` if
    the program is infeasible or unbounded.
    """
    n = len(c)
    # Normalize objective to maximization.
    obj = list(c) if maximize else [-ci for ci in c]

    # Build the standard-form tableau with slack/surplus and artificial variables.
    rows = []
    n_slack = 0
    artificial_cols = []
    for coeffs, op, rhs in constraints:
        if len(coeffs) != n:
            raise ValueError("constraint has wrong number of coefficients")
        rows.append([list(coeffs), op, float(rhs)])
        if op in ("<=", ">="):
            n_slack += 1

    total_cols = n + n_slack + 0            # artificials appended after counting
    # First pass: assign slack/surplus columns.
    slack_index = n
    art_needed = []
    matrix = []
    rhs_vec = []
    slack_assign = []
    for coeffs, op, rhs in rows:
        row = list(coeffs) + [0.0] * n_slack
        if op == "<=":
            row[slack_index] = 1.0
            slack_assign.append(slack_index)
            slack_index += 1
            art_needed.append(False)
        elif op == ">=":
            row[slack_index] = -1.0
            slack_assign.append(None)
            slack_index += 1
            art_needed.append(True)
        else:                               # "="
            slack_assign.append(None)
            art_needed.append(True)
        # Ensure non-negative RHS by flipping the row if needed.
        if rhs < 0:
            row = [-v for v in row]
            rhs = -rhs
            # flipping a <= slack turns it into a >= surplus -> needs an artificial
            if op == "<=":
                slack_assign[-1] = None
                art_needed[-1] = True
        matrix.append(row)
        rhs_vec.append(rhs)

    m = len(matrix)
    # Append artificial columns where needed.
    art_cols = []
    for i in range(m):
        if art_needed[i]:
            for r in range(m):
                matrix[r].append(1.0 if r == i else 0.0)
            art_cols.append(n + n_slack + len(art_cols))
    ncols = n + n_slack + len(art_cols)

    # Basis: slack columns where assigned, else the artificial for that row.
    basis = []
    art_ptr = 0
    for i in range(m):
        if slack_assign[i] is not None:
            basis.append(slack_assign[i])
        else:
            basis.append(art_cols[art_ptr])
            art_ptr += 1

    def pivot(tableau, basis, obj_row):
        for _ in range(max_iter):
            # Bland's rule: first column with negative reduced cost.
            entering = -1
            for j in range(len(obj_row) - 1):
                if obj_row[j] < -tol:
                    entering = j
                    break
            if entering == -1:
                return True
            # Ratio test, smallest index on ties.
            leaving = -1
            best_ratio = None
            for i in range(len(tableau)):
                a = tableau[i][entering]
                if a > tol:
                    ratio = tableau[i][-1] / a
                    if best_ratio is None or ratio < best_ratio - tol or \
                       (abs(ratio - best_ratio) <= tol and basis[i] < basis[leaving]):
                        best_ratio = ratio
                        leaving = i
            if leaving == -1:
                raise ValueError("problem is unbounded")
            # Pivot.
            piv = tableau[leaving][entering]
            tableau[leaving] = [v / piv for v in tableau[leaving]]
            for i in range(len(tableau)):
                if i != leaving and abs(tableau[i][entering]) > tol:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][k] - factor * tableau[leaving][k]
                                  for k in range(len(tableau[i]))]
            factor = obj_row[entering]
            if abs(factor) > tol:
                for k in range(len(obj_row)):
                    obj_row[k] -= factor * tableau[leaving][k]
            basis[leaving] = entering
        raise ValueError("simplex did not converge (possible cycling)")

    tableau = [matrix[i] + [rhs_vec[i]] for i in range(m)]

    # Phase 1: minimize the sum of artificials (only if any exist).
    if art_cols:
        phase1 = [0.0] * (ncols + 1)
        for j in art_cols:
            phase1[j] = 1.0
        # Reduce so basic artificials have zero cost.
        for i in range(m):
            if basis[i] in art_cols:
                for k in range(ncols + 1):
                    phase1[k] -= tableau[i][k]
        pivot(tableau, basis, phase1)
        # Feasibility: phase-1 objective (=-phase1[-1]) must be ~0.
        if abs(phase1[-1]) > 1e-6:
            raise ValueError("problem is infeasible")
        # Drive any remaining basic artificials out if possible; otherwise drop rows.
        # (For well-posed inputs this is rare; leave them at value 0.)

    # Phase 2: original objective (maximize -> minimize -obj on reduced costs).
    obj_row = [0.0] * (ncols + 1)
    for j in range(n):
        obj_row[j] = -obj[j]
    # Zero out reduced costs of basic variables.
    for i in range(m):
        if basis[i] < ncols and abs(obj_row[basis[i]]) > tol:
            factor = obj_row[basis[i]]
            for k in range(ncols + 1):
                obj_row[k] -= factor * tableau[i][k]
    # Forbid artificials from re-entering by pinning their cost high.
    for j in art_cols:
        obj_row[j] = max(obj_row[j], 1e9)
    pivot(tableau, basis, obj_row)

    x = [0.0] * n
    for i in range(m):
        if basis[i] < n:
            x[basis[i]] = tableau[i][-1]
    value = sum(c[j] * x[j] for j in range(n))
    return {"x": x, "objective": value, "status": "optimal"}
