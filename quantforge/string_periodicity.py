"""Linear-time string primitives: Z-function, prefix function, Manacher, periodicity.

Three classic ``O(n)`` string algorithms and what they unlock:

  * ``z_function`` -- ``z[i]`` is the length of the longest substring starting at ``i`` that
    matches a prefix of the string.
  * ``prefix_function`` -- ``pi[i]`` is the length of the longest proper prefix of
    ``s[:i+1]`` that is also a suffix (the KMP "failure" array); it yields the smallest
    period and the set of borders.
  * ``manacher_longest_palindrome`` -- the longest palindromic substring, found by
    Manacher's algorithm in linear time.

Pure standard library.
"""


def z_function(s):
    """Return the Z-array: ``z[i]`` = length of the longest common prefix of ``s`` and ``s[i:]``.

    ``z[0]`` is conventionally ``0``. Computed in ``O(n)`` by maintaining the rightmost
    matching segment ``[l, r]`` seen so far.
    """
    n = len(s)
    z = [0] * n
    l = r = 0
    for i in range(1, n):
        if i < r:
            z[i] = min(r - i, z[i - l])
        while i + z[i] < n and s[z[i]] == s[i + z[i]]:
            z[i] += 1
        if i + z[i] > r:
            l, r = i, i + z[i]
    return z


def prefix_function(s):
    """Return the KMP prefix function ``pi``: the longest proper prefix = suffix length."""
    n = len(s)
    pi = [0] * n
    for i in range(1, n):
        j = pi[i - 1]
        while j > 0 and s[i] != s[j]:
            j = pi[j - 1]
        if s[i] == s[j]:
            j += 1
        pi[i] = j
    return pi


def smallest_period(s):
    """Return the length of the smallest period ``p`` such that ``s`` repeats ``s[:p]``.

    A period ``p`` means ``s[i] == s[i - p]`` for all ``i >= p`` (the last block may be
    partial). For a string that is a whole number of copies this is the repeating unit; a
    string with no shorter period returns its own length. Empty string returns ``0``.
    """
    n = len(s)
    if n == 0:
        return 0
    pi = prefix_function(s)
    return n - pi[n - 1]


def is_periodic(s):
    """True if ``s`` is a whole number (>= 2) of copies of a shorter block."""
    n = len(s)
    if n == 0:
        return False
    p = smallest_period(s)
    return p < n and n % p == 0


def borders(s):
    """Return the lengths of all borders of ``s`` (proper prefix == suffix), ascending."""
    n = len(s)
    if n == 0:
        return []
    pi = prefix_function(s)
    result = []
    k = pi[n - 1]
    while k > 0:
        result.append(k)
        k = pi[k - 1]
    result.reverse()
    return result


def count_occurrences(text, pattern):
    """Return the start indices of every occurrence of ``pattern`` in ``text`` (Z-function)."""
    if pattern == "":
        return list(range(len(text) + 1))
    if len(pattern) > len(text):
        return []
    sep = "\x00"
    combined = pattern + sep + text
    z = z_function(combined)
    m = len(pattern)
    out = []
    for i in range(m + 1, len(combined)):
        if z[i] >= m:
            out.append(i - m - 1)
    return out


def manacher_longest_palindrome(s):
    """Return a longest palindromic substring of ``s`` (Manacher's algorithm, ``O(n)``).

    Ties resolve to the earliest-starting longest palindrome. The empty string returns
    ``""``.
    """
    if s == "":
        return ""
    # transform to handle even/odd uniformly: ^ # a # b # ... $
    t = "^#" + "#".join(s) + "#$"
    n = len(t)
    p = [0] * n
    center = right = 0
    for i in range(1, n - 1):
        if i < right:
            p[i] = min(right - i, p[2 * center - i])
        while t[i + p[i] + 1] == t[i - p[i] - 1]:
            p[i] += 1
        if i + p[i] > right:
            center, right = i, i + p[i]
    # find the max radius; earliest on ties
    best_len = 0
    best_center = 0
    for i in range(1, n - 1):
        if p[i] > best_len:
            best_len = p[i]
            best_center = i
    start = (best_center - best_len) // 2
    return s[start:start + best_len]
