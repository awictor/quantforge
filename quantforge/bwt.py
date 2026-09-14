"""Burrows-Wheeler transform and move-to-front coding (bzip2 building blocks).

The Burrows-Wheeler transform reversibly permutes a string so that runs of the same
character cluster together, making it far more compressible -- it is the heart of bzip2.
Move-to-front coding then turns those clusters into runs of small integers (mostly
zeros), ideal for a following run-length + entropy coder. Both invert exactly. A sentinel
index (rather than an appended terminator) makes the BWT self-inverting on arbitrary
input. Pure standard library.
"""


def bwt_transform(data):
    """Burrows-Wheeler transform of a sequence: ``(transformed, primary_index)``.

    Builds the sorted matrix of all rotations and returns the last column plus the index
    of the original string among the sorted rotations (needed to invert). ``data`` is a
    string or list; the transformed output has the same type. Empty input returns
    ``("" or [], 0)``.
    """
    n = len(data)
    if n == 0:
        return (data[:], 0)
    # Rotations as index offsets; sort by the rotated sequence.
    rotations = sorted(range(n), key=lambda i: data[i:] + data[:i])
    is_str = isinstance(data, str)
    last = [data[(i + n - 1) % n] for i in rotations]
    primary = rotations.index(0)
    return ("".join(last) if is_str else last, primary)


def bwt_inverse(transformed, primary_index):
    """Invert the Burrows-Wheeler transform back to the original sequence.

    ``transformed`` is the last column and ``primary_index`` the row of the original
    string (both from :func:`bwt_transform`). Reconstructs via the standard LF-mapping.
    Returns the same type as ``transformed``.
    """
    n = len(transformed)
    if n == 0:
        return transformed[:]
    is_str = isinstance(transformed, str)
    last = list(transformed)
    # Stable sort of (symbol, original position) is the LF-mapping (left-shift permutation).
    lshift = sorted(range(n), key=lambda i: (last[i], i))
    # Walk the mapping from the primary row, stepping before reading.
    out = []
    row = primary_index
    for _ in range(n):
        row = lshift[row]
        out.append(last[row])
    return "".join(out) if is_str else out


def move_to_front_encode(data, alphabet=None):
    """Move-to-front encode a sequence into a list of integer ranks.

    Each symbol is replaced by its current index in a running alphabet list, then moved
    to the front. Clustered inputs (like BWT output) produce many zeros. ``alphabet`` is
    the initial ordered symbol list; by default the sorted set of symbols in ``data``.
    Returns ``(codes, alphabet)`` -- the alphabet is needed to decode.
    """
    if alphabet is None:
        alphabet = sorted(set(data))
    table = list(alphabet)
    codes = []
    for s in data:
        idx = table.index(s)
        codes.append(idx)
        table.pop(idx)
        table.insert(0, s)
    return codes, list(alphabet)


def move_to_front_decode(codes, alphabet):
    """Invert move-to-front coding given the ``codes`` and the initial ``alphabet``.

    Returns the reconstructed sequence: a ``str`` if the alphabet symbols are all
    single characters, else a list.
    """
    table = list(alphabet)
    out = []
    for idx in codes:
        s = table[idx]
        out.append(s)
        table.pop(idx)
        table.insert(0, s)
    if all(isinstance(s, str) and len(s) == 1 for s in alphabet):
        return "".join(out)
    return out
