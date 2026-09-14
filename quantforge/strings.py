"""String algorithms: edit distance, common subsequence/substring, pattern search.

Classic sequence-comparison and matching routines: Levenshtein edit distance and its
Hamming special case, the longest common subsequence and substring (dynamic
programming), and Knuth-Morris-Pratt substring search in linear time. These drive fuzzy
matching, diffing, and text search. Work on any sequences (strings or lists). Pure
standard library.
"""


def levenshtein(a, b):
    """Levenshtein edit distance: minimum single-character insert/delete/substitutions.

    The number of one-character edits to turn ``a`` into ``b``. ``O(len(a)*len(b))`` time
    with two rolling rows. Symmetric; zero iff the sequences are equal.
    """
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cost = 0 if ca == cb else 1
            cur.append(min(prev[j] + 1,          # deletion
                           cur[j - 1] + 1,        # insertion
                           prev[j - 1] + cost))   # substitution
        prev = cur
    return prev[len(b)]


def hamming_distance(a, b):
    """Hamming distance: positions at which equal-length sequences differ.

    Raises ``ValueError`` if the lengths differ (Hamming distance is undefined then).
    """
    if len(a) != len(b):
        raise ValueError("hamming_distance requires equal-length sequences")
    return sum(1 for x, y in zip(a, b) if x != y)


def longest_common_subsequence(a, b):
    """Longest common subsequence (not necessarily contiguous) as a list/str.

    Returns a longest sequence appearing in both ``a`` and ``b`` in order. Returns the
    same type as ``a`` when ``a`` is a ``str``. ``O(len(a)*len(b))``.
    """
    m, n = len(a), len(b)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if a[i - 1] == b[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])
    # Backtrack.
    out = []
    i, j = m, n
    while i > 0 and j > 0:
        if a[i - 1] == b[j - 1]:
            out.append(a[i - 1])
            i -= 1
            j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            i -= 1
        else:
            j -= 1
    out.reverse()
    return "".join(out) if isinstance(a, str) else out


def longest_common_substring(a, b):
    """Longest contiguous substring common to ``a`` and ``b``.

    Returns the substring (same type as ``a`` when ``a`` is a ``str``); an empty result
    if there is no common character. ``O(len(a)*len(b))``.
    """
    m, n = len(a), len(b)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    best_len = 0
    best_end = 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if a[i - 1] == b[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
                if dp[i][j] > best_len:
                    best_len = dp[i][j]
                    best_end = i
    result = a[best_end - best_len:best_end]
    return result


def kmp_search(text, pattern):
    """All start indices where ``pattern`` occurs in ``text`` (Knuth-Morris-Pratt).

    Linear-time ``O(len(text)+len(pattern))`` search using the failure function, so it
    never re-examines text characters. An empty pattern matches at every position
    ``0..len(text)``. Returns a list of indices.
    """
    if not pattern:
        return list(range(len(text) + 1))
    # Build the failure (longest proper prefix-suffix) table.
    lps = [0] * len(pattern)
    k = 0
    for i in range(1, len(pattern)):
        while k > 0 and pattern[i] != pattern[k]:
            k = lps[k - 1]
        if pattern[i] == pattern[k]:
            k += 1
        lps[i] = k
    matches = []
    k = 0
    for i in range(len(text)):
        while k > 0 and text[i] != pattern[k]:
            k = lps[k - 1]
        if text[i] == pattern[k]:
            k += 1
        if k == len(pattern):
            matches.append(i - k + 1)
            k = lps[k - 1]
    return matches
