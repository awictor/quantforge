"""Lossless compression: Huffman coding and run-length encoding."""

import math
import random
from collections import Counter

import pytest

from quantforge import (
    huffman_codebook,
    huffman_encode,
    huffman_decode,
    run_length_encode,
    run_length_decode,
)


def test_huffman_roundtrip():
    rng = random.Random(1)
    for _ in range(300):
        data = [rng.choice("abcde") for _ in range(rng.randint(1, 50))]
        bits, cb = huffman_encode(data)
        assert huffman_decode(bits, cb) == data


def test_huffman_prefix_free():
    rng = random.Random(2)
    for _ in range(200):
        data = [rng.choice("abcdefgh") for _ in range(rng.randint(1, 40))]
        codes = list(huffman_codebook(data).values())
        for i in range(len(codes)):
            for j in range(len(codes)):
                if i != j:
                    assert not codes[i].startswith(codes[j])


def test_huffman_near_entropy_and_beats_fixed():
    data = "mississippi" * 10
    bits, _ = huffman_encode(data)
    n = len(data)
    freq = Counter(data)
    H = -sum((f / n) * math.log2(f / n) for f in freq.values())
    avg = len(bits) / n
    assert H <= avg < H + 1                      # Huffman within [H, H+1)
    assert avg <= math.ceil(math.log2(len(freq)))  # beats fixed-length


def test_huffman_single_symbol():
    bits, cb = huffman_encode("aaaa")
    assert cb == {"a": "0"}
    assert "".join(huffman_decode(bits, cb)) == "aaaa"


def test_huffman_frequent_symbol_shortest():
    cb = huffman_codebook("aaaabbbccd")           # a:4 b:3 c:2 d:1
    assert len(cb["a"]) <= len(cb["d"])


def test_rle_roundtrip():
    rng = random.Random(3)
    for _ in range(300):
        data = [rng.choice("ab") for _ in range(rng.randint(0, 40))]
        assert run_length_decode(run_length_encode(data)) == data
    assert run_length_encode("aaabbc") == [("a", 3), ("b", 2), ("c", 1)]


def test_validation():
    with pytest.raises(ValueError):
        huffman_encode("")
    with pytest.raises(ValueError):
        run_length_decode([("a", 0)])
    with pytest.raises(ValueError):
        huffman_decode("111", {"a": "0"})
