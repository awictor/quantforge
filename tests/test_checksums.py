"""Checksums and non-cryptographic hashes: CRC-32, Adler-32, FNV-1a."""

import random
import zlib

from quantforge import crc32, adler32, fnv1a_32


def test_crc32_matches_zlib():
    rng = random.Random(1)
    for _ in range(1000):
        b = bytes(rng.randint(0, 255) for _ in range(rng.randint(0, 50)))
        assert crc32(b) == zlib.crc32(b)
    assert crc32("hello") == zlib.crc32(b"hello")
    assert crc32(b"") == 0


def test_adler32_matches_zlib():
    rng = random.Random(2)
    for _ in range(1000):
        b = bytes(rng.randint(0, 255) for _ in range(rng.randint(0, 50)))
        assert adler32(b) == zlib.adler32(b)
    assert adler32("Wikipedia") == 300286872


def test_fnv1a_known_vectors():
    assert fnv1a_32("") == 2166136261
    assert fnv1a_32("a") == 3826002220
    assert fnv1a_32("foobar") == 3214735720


def test_str_bytes_equivalence():
    assert crc32("hello") == crc32(b"hello")
    assert adler32("hello") == adler32(b"hello")
    assert fnv1a_32("test") == fnv1a_32(b"test")


def test_outputs_are_32bit():
    for fn in (crc32, adler32, fnv1a_32):
        v = fn("the quick brown fox")
        assert 0 <= v <= 0xFFFFFFFF
