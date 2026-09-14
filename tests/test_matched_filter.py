"""Matched filtering and peak detection."""

import math
import random

import pytest

from quantforge import (
    matched_filter,
    normalized_matched_filter,
    find_peaks,
    detect_template,
)

TEMPLATE = [1.0, 2.0, 3.0, 2.0, 1.0]


def test_matched_filter_peaks_at_location():
    rng = random.Random(5)
    N, loc = 200, 80
    x = [rng.gauss(0.0, 0.5) for _ in range(N)]
    for j in range(len(TEMPLATE)):
        x[loc + j] += TEMPLATE[j] * 3.0
    mf = matched_filter(x, TEMPLATE)
    assert max(range(len(mf)), key=lambda i: mf[i]) == loc
    assert len(mf) == N - len(TEMPLATE) + 1


def test_normalized_is_one_at_scaled_match():
    clean = [0.0] * 50
    for j in range(len(TEMPLATE)):
        clean[20 + j] = TEMPLATE[j] * 7.0
    nm = normalized_matched_filter(clean, TEMPLATE)
    assert abs(max(nm) - 1.0) < 1e-9
    assert max(range(len(nm)), key=lambda i: nm[i]) == 20


def test_find_peaks_basic_and_height():
    sig = [0, 3, 0, 5, 0, 1, 0, 9, 0]
    assert find_peaks(sig) == [1, 3, 5, 7]
    assert find_peaks(sig, height=4) == [3, 7]


def test_find_peaks_distance_keeps_taller():
    close = [0, 5, 0, 6, 0]        # peaks at index 1 and 3
    assert find_peaks(close, distance=2) == [1, 3]
    tight = [0, 5, 6, 0]           # not both strict local maxima
    assert find_peaks(tight, distance=3) == [2]


def test_detect_template_finds_all_copies():
    # a longer distinctive template is unambiguous even in noise
    template = [math.sin(math.pi * k / 16) * math.cos(2 * math.pi * 3 * k / 16)
                for k in range(16)]
    rng = random.Random(11)
    N = 400
    y = [rng.gauss(0.0, 0.2) for _ in range(N)]
    locs = [50, 160, 300]
    for L in locs:
        for j in range(len(template)):
            y[L + j] += template[j] * 3.0
    assert detect_template(y, template, threshold=0.9) == locs


def test_validation():
    with pytest.raises(ValueError):
        matched_filter([], [1.0])
    with pytest.raises(ValueError):
        matched_filter([1.0, 2.0], [1.0, 2.0, 3.0])     # template longer than signal
    with pytest.raises(ValueError):
        normalized_matched_filter([1.0, 2.0, 3.0], [5.0, 5.0])   # zero-variance template
    with pytest.raises(ValueError):
        find_peaks([1.0, 2.0, 1.0], distance=0)
