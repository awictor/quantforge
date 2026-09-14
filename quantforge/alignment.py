"""Sequence alignment: Damerau-Levenshtein distance and Needleman-Wunsch.

Extends plain edit distance with a *transposition* operation (swapping two adjacent
characters counts as one edit, not two) -- the true (unrestricted) Damerau-Levenshtein
distance, which matches how typos actually occur. Needleman-Wunsch is the classic global
alignment: it maximizes a match/mismatch/gap score and reconstructs the aligned strings.
Pure standard library.
"""


def damerau_levenshtein(a, b):
    """True (unrestricted) Damerau-Levenshtein distance between two sequences.

    Counts insertions, deletions, substitutions, and transpositions of adjacent elements,
    each as a single edit -- so ``"ca" -> "ac"`` is distance 1, not 2. Uses the full
    dynamic-programming table with a last-seen index (unlike the restricted OSA variant).
    """
    la, lb = len(a), len(b)
    max_dist = la + lb
    # Alphabet of symbols seen, mapped to a "last row" tracker.
    da = {}
    # DP table padded by one row/col of sentinels holding max_dist.
    d = [[0] * (lb + 2) for _ in range(la + 2)]
    d[0][0] = max_dist
    for i in range(la + 1):
        d[i + 1][0] = max_dist
        d[i + 1][1] = i
    for j in range(lb + 1):
        d[0][j + 1] = max_dist
        d[1][j + 1] = j
    for i in range(1, la + 1):
        db = 0
        for j in range(1, lb + 1):
            k = da.get(b[j - 1], 0)
            ell = db
            if a[i - 1] == b[j - 1]:
                cost = 0
                db = j
            else:
                cost = 1
            d[i + 1][j + 1] = min(
                d[i][j] + cost,                              # substitution
                d[i + 1][j] + 1,                             # insertion
                d[i][j + 1] + 1,                             # deletion
                d[k][ell] + (i - k - 1) + 1 + (j - ell - 1)  # transposition
            )
        da[a[i - 1]] = i
    return d[la + 1][lb + 1]


def needleman_wunsch(a, b, match=1, mismatch=-1, gap=-1):
    """Needleman-Wunsch global alignment: ``(score, aligned_a, aligned_b)``.

    Maximizes the total score with ``match``/``mismatch`` for aligned pairs and ``gap``
    per inserted gap ``'-'``. Returns the optimal score and the two gapped strings (as
    ``str`` when the inputs are strings, else lists). Standard ``O(len(a)*len(b))`` DP with
    traceback.
    """
    la, lb = len(a), len(b)
    # Score matrix.
    f = [[0] * (lb + 1) for _ in range(la + 1)]
    for i in range(1, la + 1):
        f[i][0] = i * gap
    for j in range(1, lb + 1):
        f[0][j] = j * gap
    for i in range(1, la + 1):
        for j in range(1, lb + 1):
            s = match if a[i - 1] == b[j - 1] else mismatch
            f[i][j] = max(f[i - 1][j - 1] + s, f[i - 1][j] + gap, f[i][j - 1] + gap)
    # Traceback.
    ai, bj = [], []
    i, j = la, lb
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = match if a[i - 1] == b[j - 1] else mismatch
            if f[i][j] == f[i - 1][j - 1] + s:
                ai.append(a[i - 1]); bj.append(b[j - 1]); i -= 1; j -= 1; continue
        if i > 0 and f[i][j] == f[i - 1][j] + gap:
            ai.append(a[i - 1]); bj.append("-"); i -= 1
        else:
            ai.append("-"); bj.append(b[j - 1]); j -= 1
    ai.reverse(); bj.reverse()
    is_str = isinstance(a, str) and isinstance(b, str)
    if is_str:
        return f[la][lb], "".join(ai), "".join(bj)
    return f[la][lb], ai, bj
