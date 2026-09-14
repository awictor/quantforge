"""Great-circle geodesy on a spherical Earth."""

import random

import pytest

from quantforge import (
    haversine_distance,
    initial_bearing,
    destination_point,
    cross_track_distance,
)


def test_known_city_distances():
    assert abs(haversine_distance(51.5074, -0.1278, 48.8566, 2.3522) - 343.6) < 2
    assert abs(haversine_distance(40.7128, -74.0060, 34.0522, -118.2437) - 3936) < 5


def test_distance_symmetric_and_zero():
    rng = random.Random(1)
    for _ in range(500):
        a, b, c, d = (rng.uniform(-80, 80), rng.uniform(-180, 180),
                      rng.uniform(-80, 80), rng.uniform(-180, 180))
        assert abs(haversine_distance(a, b, c, d) - haversine_distance(c, d, a, b)) < 1e-9
    assert haversine_distance(45, 45, 45, 45) == 0.0


def test_destination_roundtrip():
    rng = random.Random(2)
    for _ in range(500):
        lat, lon = rng.uniform(-70, 70), rng.uniform(-170, 170)
        brng, dist = rng.uniform(0, 360), rng.uniform(1, 5000)
        lat2, lon2 = destination_point(lat, lon, brng, dist)
        assert abs(haversine_distance(lat, lon, lat2, lon2) - dist) < 1e-6
        b2 = initial_bearing(lat, lon, lat2, lon2)
        assert abs((b2 - brng + 180) % 360 - 180) < 1e-4


def test_cardinal_bearings():
    assert abs(initial_bearing(0, 0, 10, 0)) < 1e-9        # due north
    assert abs(initial_bearing(0, 0, 0, 10) - 90) < 1e-9   # due east


def test_cross_track():
    assert abs(cross_track_distance(5, 0, 0, 0, 10, 0)) < 1e-6      # on the meridian
    off = cross_track_distance(5, 1, 0, 0, 10, 0)                   # 1 deg east
    assert abs(abs(off) - 111) < 3                                  # ~111 km per degree
