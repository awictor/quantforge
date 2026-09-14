"""LZW dictionary compression and integer delta coding.

LZW builds its dictionary of repeated substrings on the fly -- no separate codebook is
transmitted, and it adapts to the data as it reads. Delta coding replaces a sequence of
numbers by their successive differences, shrinking slowly-varying series so a following
entropy coder sees small values. Both invert losslessly. Pure standard library.
"""


def lzw_compress(data):
    """LZW-compress a string into a list of integer codes.

    Starts with a dictionary of the distinct single characters (assigned codes in sorted
    order) and grows it with each new substring seen, emitting the code of the longest
    known prefix. Returns ``(codes, alphabet)`` -- the sorted initial alphabet is needed
    to seed the decoder. Empty input yields ``([], [])``.
    """
    if not data:
        return [], []
    alphabet = sorted(set(data))
    table = {ch: i for i, ch in enumerate(alphabet)}
    next_code = len(alphabet)
    codes = []
    w = ""
    for ch in data:
        wc = w + ch
        if wc in table:
            w = wc
        else:
            codes.append(table[w])
            table[wc] = next_code
            next_code += 1
            w = ch
    if w:
        codes.append(table[w])
    return codes, alphabet


def lzw_decompress(codes, alphabet):
    """Decompress LZW ``codes`` given the initial ``alphabet`` back to the string.

    Rebuilds the same dictionary the compressor grew, handling the special case where a
    code refers to an entry being defined this step. Inverts :func:`lzw_compress`
    exactly.
    """
    if not codes:
        return ""
    table = {i: ch for i, ch in enumerate(alphabet)}
    next_code = len(alphabet)
    if codes[0] not in table:
        raise ValueError("invalid LZW code stream")
    w = table[codes[0]]
    out = [w]
    for code in codes[1:]:
        if code in table:
            entry = table[code]
        elif code == next_code:
            entry = w + w[0]          # the not-yet-defined case
        else:
            raise ValueError("invalid LZW code stream")
        out.append(entry)
        table[next_code] = w + entry[0]
        next_code += 1
        w = entry
    return "".join(out)


def delta_encode(data):
    """Delta-encode a numeric sequence: first value, then successive differences.

    ``[a, b, c, ...] -> [a, b-a, c-b, ...]``. Slowly-varying data becomes small numbers
    that compress better. Empty input yields an empty list.
    """
    out = []
    prev = 0
    for i, v in enumerate(data):
        out.append(v - prev if i else v)
        prev = v
    return out


def delta_decode(deltas):
    """Invert :func:`delta_encode`: cumulative sum of the deltas."""
    out = []
    total = 0
    for d in deltas:
        total += d
        out.append(total)
    return out
