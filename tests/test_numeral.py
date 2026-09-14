"""Numeral systems: base conversion and Roman numerals."""

import random

import pytest

from quantforge import to_base, from_base, to_roman, from_roman


def test_base_roundtrip_vs_int():
    rng = random.Random(1)
    for _ in range(5000):
        n = rng.randint(0, 10 ** 9)
        base = rng.randint(2, 36)
        s = to_base(n, base)
        assert from_base(s, base) == n
        assert int(s, base) == n


def test_base_known():
    assert to_base(255, 16) == "ff"
    assert to_base(10, 2) == "1010"
    assert to_base(0, 7) == "0"
    assert from_base("FF", 16) == 255          # case-insensitive


def test_roman_known():
    known = {1: "I", 4: "IV", 9: "IX", 40: "XL", 90: "XC", 400: "CD", 900: "CM",
             1994: "MCMXCIV", 2024: "MMXXIV", 3999: "MMMCMXCIX"}
    for n, r in known.items():
        assert to_roman(n) == r
        assert from_roman(r) == n


def test_roman_roundtrip():
    for n in range(1, 4000):
        assert from_roman(to_roman(n)) == n


def test_validation():
    with pytest.raises(ValueError):
        to_base(-1, 10)
    with pytest.raises(ValueError):
        to_base(5, 40)
    with pytest.raises(ValueError):
        from_base("2", 2)
    with pytest.raises(ValueError):
        to_roman(0)
    with pytest.raises(ValueError):
        to_roman(4000)
    with pytest.raises(ValueError):
        from_roman("ABC")
