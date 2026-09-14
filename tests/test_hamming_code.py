"""Error-correcting codes: Hamming(7,4) and the Luhn checksum."""

import random

import pytest

from quantforge import (
    hamming74_encode,
    hamming74_decode,
    luhn_checksum,
    luhn_check_digit,
)


def test_hamming_clean_roundtrip():
    for x in range(16):
        bits = [(x >> i) & 1 for i in range(4)]
        data, err = hamming74_decode(hamming74_encode(bits))
        assert data == bits and err == 0


def test_hamming_corrects_single_bit_errors():
    for x in range(16):
        bits = [(x >> i) & 1 for i in range(4)]
        code = hamming74_encode(bits)
        for pos in range(7):
            corrupt = list(code)
            corrupt[pos] ^= 1
            data, err = hamming74_decode(corrupt)
            assert data == bits            # corrected
            assert err == pos + 1          # located the flip


def test_luhn_known_numbers():
    assert luhn_checksum("79927398713") == 0
    assert luhn_checksum("4532015112830366") == 0     # valid test Visa
    assert luhn_checksum("79927398710") != 0


def test_luhn_check_digit_consistency():
    rng = random.Random(1)
    for _ in range(1000):
        payload = [rng.randint(0, 9) for _ in range(rng.randint(1, 15))]
        cd = luhn_check_digit(payload)
        assert luhn_checksum(payload + [cd]) == 0
    assert luhn_check_digit("7992739871") == 3


def test_validation():
    with pytest.raises(ValueError):
        hamming74_encode([0, 1, 2, 0])
    with pytest.raises(ValueError):
        hamming74_decode([0, 0, 0])
