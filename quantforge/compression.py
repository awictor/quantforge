"""Lossless compression primitives: Huffman coding and run-length encoding.

Two foundational lossless codes. Huffman builds the optimal prefix code for a symbol
distribution (shortest expected code length of any prefix code), returning the bitstring
and codebook so it round-trips exactly. Run-length encoding collapses runs of a repeated
symbol into ``(symbol, count)`` pairs -- ideal for sparse or blocky data. Both invert
losslessly. Pure standard library.
"""

import heapq
from collections import Counter


def huffman_codebook(data):
    """Optimal Huffman prefix-code table ``{symbol: bitstring}`` for a sequence.

    Builds the code by repeatedly merging the two least-frequent nodes. A single distinct
    symbol maps to ``"0"`` (a one-bit code). Raises on empty input. The codes are
    prefix-free, so no code is a prefix of another.
    """
    if not data:
        raise ValueError("input must be non-empty")
    freq = Counter(data)
    if len(freq) == 1:
        return {next(iter(freq)): "0"}
    # Min-heap of (frequency, tie-breaker, subtree). Subtree is a leaf symbol or a pair.
    heap = [(f, i, sym) for i, (sym, f) in enumerate(sorted(freq.items(), key=lambda kv: (kv[1], repr(kv[0]))))]
    heapq.heapify(heap)
    counter = len(heap)
    while len(heap) > 1:
        f1, _, left = heapq.heappop(heap)
        f2, _, right = heapq.heappop(heap)
        heapq.heappush(heap, (f1 + f2, counter, (left, right)))
        counter += 1
    root = heap[0][2]
    codebook = {}

    def walk(node, prefix):
        if isinstance(node, tuple):
            walk(node[0], prefix + "0")
            walk(node[1], prefix + "1")
        else:
            codebook[node] = prefix

    walk(root, "")
    return codebook


def huffman_encode(data):
    """Huffman-encode a sequence: return ``(bitstring, codebook)``.

    ``bitstring`` is a ``str`` of ``'0'``/``'1'``; ``codebook`` maps each symbol to its
    code (needed to decode). Expected length is minimal among prefix codes.
    """
    codebook = huffman_codebook(data)
    bits = "".join(codebook[s] for s in data)
    return bits, codebook


def huffman_decode(bits, codebook):
    """Decode a Huffman ``bits`` string given its ``codebook`` back to the symbol list.

    Inverts :func:`huffman_encode` exactly. Raises if the bitstring is not a valid
    concatenation of codes.
    """
    inverse = {code: sym for sym, code in codebook.items()}
    out = []
    cur = ""
    for b in bits:
        cur += b
        if cur in inverse:
            out.append(inverse[cur])
            cur = ""
    if cur:
        raise ValueError("bitstring is not a valid Huffman encoding")
    return out


def run_length_encode(data):
    """Run-length encode a sequence into a list of ``(symbol, count)`` pairs.

    Consecutive equal symbols collapse into one pair; ``count >= 1``. Empty input yields
    an empty list.
    """
    out = []
    for s in data:
        if out and out[-1][0] == s:
            out[-1] = (s, out[-1][1] + 1)
        else:
            out.append((s, 1))
    return out


def run_length_decode(pairs):
    """Expand ``(symbol, count)`` pairs back into the original sequence (a list).

    Inverts :func:`run_length_encode`. Raises on a non-positive count.
    """
    out = []
    for sym, count in pairs:
        if count < 1:
            raise ValueError("count must be >= 1")
        out.extend([sym] * count)
    return out
