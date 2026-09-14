"""LZW dictionary compression and integer delta coding."""

import random

import pytest

from quantforge import lzw_compress, lzw_decompress, delta_encode, delta_decode


def test_lzw_roundtrip_known():
    for s in ["TOBEORNOTTOBEORTOBEORNOT", "", "a", "aaaaaa",
              "ababababab", "abcabcabc"]:
        codes, alpha = lzw_compress(s)
        assert lzw_decompress(codes, alpha) == s


def test_lzw_roundtrip_random():
    rng = random.Random(1)
    for _ in range(500):
        s = "".join(rng.choice("ab") for _ in range(rng.randint(0, 60)))
        codes, alpha = lzw_compress(s)
        assert lzw_decompress(codes, alpha) == s


def test_lzw_compresses_repetitive():
    s = "abcabcabcabcabcabcabcabc"
    codes, _ = lzw_compress(s)
    assert len(codes) < len(s)


def test_delta_roundtrip():
    rng = random.Random(2)
    for _ in range(500):
        data = [rng.randint(-100, 100) for _ in range(rng.randint(0, 30))]
        assert delta_decode(delta_encode(data)) == data


def test_delta_known():
    assert delta_encode([10, 12, 15, 14]) == [10, 2, 3, -1]
    assert delta_encode([100, 101, 102, 103, 104]) == [100, 1, 1, 1, 1]


def test_validation():
    with pytest.raises(ValueError):
        lzw_decompress([999], ["a"])
