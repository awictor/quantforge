"""Fuzzy string similarity: Jaro, Jaro-Winkler, Dice, Jaccard.

Similarity scores in ``[0, 1]`` for approximate string matching -- deduplication, record
linkage, spell correction, and search ranking. Jaro (and its Winkler prefix boost) is
tuned for short strings like names; the Dice and Jaccard coefficients compare character
bigram (or token) sets. All are ``1`` for identical inputs and ``0`` for wholly
dissimilar ones. Pure standard library.
"""


def jaro(a, b):
    """Jaro similarity of two strings in ``[0, 1]`` (1 = identical).

    Counts matching characters within a sliding window and penalizes transpositions,
    the classic short-string metric. Empty-vs-empty is ``1``; empty-vs-nonempty is ``0``.
    """
    if a == b:
        return 1.0
    la, lb = len(a), len(b)
    if la == 0 or lb == 0:
        return 0.0
    window = max(la, lb) // 2 - 1
    if window < 0:
        window = 0
    a_match = [False] * la
    b_match = [False] * lb
    matches = 0
    for i in range(la):
        lo = max(0, i - window)
        hi = min(i + window + 1, lb)
        for j in range(lo, hi):
            if not b_match[j] and a[i] == b[j]:
                a_match[i] = True
                b_match[j] = True
                matches += 1
                break
    if matches == 0:
        return 0.0
    # Count transpositions.
    transpositions = 0
    k = 0
    for i in range(la):
        if a_match[i]:
            while not b_match[k]:
                k += 1
            if a[i] != b[k]:
                transpositions += 1
            k += 1
    transpositions //= 2
    m = matches
    return (m / la + m / lb + (m - transpositions) / m) / 3.0


def jaro_winkler(a, b, prefix_weight=0.1, max_prefix=4):
    """Jaro-Winkler similarity: Jaro boosted for a shared prefix.

    Adds ``prefix_weight * L * (1 - jaro)`` where ``L`` is the common prefix length (up
    to ``max_prefix``), rewarding strings that agree at the start -- the standard tweak
    for names. ``prefix_weight`` must satisfy ``prefix_weight * max_prefix <= 1``.
    """
    j = jaro(a, b)
    prefix = 0
    for x, y in zip(a, b):
        if x == y and prefix < max_prefix:
            prefix += 1
        else:
            break
    return j + prefix * prefix_weight * (1.0 - j)


def _bigrams(s):
    return [s[i:i + 2] for i in range(len(s) - 1)]


def dice_coefficient(a, b):
    """Sørensen-Dice coefficient over character bigrams: ``2|A∩B| / (|A|+|B|)``.

    A multiset overlap of adjacent character pairs -- robust to word-order and small
    edits. Identical strings score ``1``; strings sharing no bigram score ``0``. Strings
    shorter than two characters fall back to exact equality.
    """
    if a == b:
        return 1.0
    if len(a) < 2 or len(b) < 2:
        return 0.0
    ba = _bigrams(a)
    bb = _bigrams(b)
    # Multiset intersection.
    from collections import Counter
    inter = sum((Counter(ba) & Counter(bb)).values())
    return 2.0 * inter / (len(ba) + len(bb))


def jaccard_similarity(a, b, tokenize=None):
    """Jaccard similarity ``|A∩B| / |A∪B|`` over token (or character) sets.

    By default splits on whitespace into a word set; pass ``tokenize`` (a callable
    returning an iterable of tokens, e.g. ``list`` for characters) to change the
    granularity. Two empty token sets score ``1``.
    """
    if tokenize is None:
        set_a, set_b = set(a.split()), set(b.split())
    else:
        set_a, set_b = set(tokenize(a)), set(tokenize(b))
    if not set_a and not set_b:
        return 1.0
    union = set_a | set_b
    if not union:
        return 1.0
    return len(set_a & set_b) / len(union)
