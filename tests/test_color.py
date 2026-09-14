"""Color-space conversions: RGB, HSV, HSL, hex."""

import colorsys
import random

import pytest

from quantforge import (
    rgb_to_hsv,
    hsv_to_rgb,
    rgb_to_hsl,
    hsl_to_rgb,
    rgb_to_hex,
    hex_to_rgb,
)


def _close(a, b, t=1e-9):
    return all(abs(a[i] - b[i]) < t for i in range(len(a)))


def test_hsv_matches_colorsys():
    rng = random.Random(1)
    for _ in range(2000):
        r, g, b = rng.random(), rng.random(), rng.random()
        h, s, v = rgb_to_hsv(r, g, b)
        ch, cs, cv = colorsys.rgb_to_hsv(r, g, b)
        assert abs(h / 360 - ch) < 1e-9 and abs(s - cs) < 1e-9 and abs(v - cv) < 1e-9


def test_hsl_matches_colorsys():
    rng = random.Random(2)
    for _ in range(2000):
        r, g, b = rng.random(), rng.random(), rng.random()
        h, s, l = rgb_to_hsl(r, g, b)
        ch, cl, cs = colorsys.rgb_to_hls(r, g, b)   # colorsys order is h, l, s
        assert abs(h / 360 - ch) < 1e-9 and abs(s - cs) < 1e-9 and abs(l - cl) < 1e-9


def test_roundtrips():
    rng = random.Random(3)
    for _ in range(2000):
        rgb = (rng.random(), rng.random(), rng.random())
        assert _close(hsv_to_rgb(*rgb_to_hsv(*rgb)), rgb)
        assert _close(hsl_to_rgb(*rgb_to_hsl(*rgb)), rgb)


def test_known_colors():
    assert _close(rgb_to_hsv(1, 0, 0), (0.0, 1.0, 1.0))
    assert _close(rgb_to_hsv(0, 1, 0), (120.0, 1.0, 1.0))
    assert _close(rgb_to_hsl(0, 0, 1), (240.0, 1.0, 0.5))


def test_hex():
    assert rgb_to_hex(1, 1, 1) == "#ffffff"
    assert rgb_to_hex(0, 0, 0) == "#000000"
    assert rgb_to_hex(1, 0, 0) == "#ff0000"
    rng = random.Random(4)
    for _ in range(500):
        code = "#%06x" % rng.randint(0, 0xFFFFFF)
        assert rgb_to_hex(*hex_to_rgb(code)) == code


def test_validation():
    with pytest.raises(ValueError):
        hex_to_rgb("#fff")
    with pytest.raises(ValueError):
        hex_to_rgb("#gggggg")
