"""Burrows-Wheeler transform and move-to-front coding."""

import random

from quantforge import (
    bwt_transform,
    bwt_inverse,
    move_to_front_encode,
    move_to_front_decode,
)


def test_bwt_roundtrip_known():
    for s in ["banana", "mississippi", "abracadabra", "", "a", "aaaa",
              "the quick brown fox"]:
        t, pi = bwt_transform(s)
        assert bwt_inverse(t, pi) == s


def test_bwt_roundtrip_random():
    rng = random.Random(1)
    for _ in range(500):
        s = "".join(rng.choice("abc") for _ in range(rng.randint(0, 30)))
        t, pi = bwt_transform(s)
        assert bwt_inverse(t, pi) == s


def test_bwt_banana():
    assert bwt_transform("banana") == ("nnbaaa", 3)


def test_mtf_roundtrip():
    rng = random.Random(2)
    for _ in range(500):
        s = "".join(rng.choice("abcd") for _ in range(rng.randint(1, 30)))
        codes, alpha = move_to_front_encode(s)
        assert move_to_front_decode(codes, alpha) == s


def test_mtf_clusters_to_zeros():
    codes, _ = move_to_front_encode("aaaabbbbcccc")
    assert codes.count(0) == 10         # only the first of b and c runs is non-zero


def test_bwt_mtf_pipeline():
    rng = random.Random(3)
    for _ in range(200):
        s = "".join(rng.choice("ab") for _ in range(rng.randint(1, 40)))
        t, pi = bwt_transform(s)
        codes, alpha = move_to_front_encode(t)
        assert bwt_inverse(move_to_front_decode(codes, alpha), pi) == s


def test_mtf_list_input():
    codes, alpha = move_to_front_encode([3, 3, 1, 1, 2])
    assert move_to_front_decode(codes, alpha) == [3, 3, 1, 1, 2]
