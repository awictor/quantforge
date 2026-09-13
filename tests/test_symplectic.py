"""Symplectic integrators: velocity Verlet and leapfrog."""

import math

import pytest

from quantforge import velocity_verlet, leapfrog


def test_harmonic_energy_bounded():
    qs, ps = velocity_verlet(lambda q: -q, 1.0, 0.0, 1.0, 0.05, 2000)
    E = [0.5 * (ps[i][0] ** 2 + qs[i][0] ** 2) for i in range(len(qs))]
    assert abs(E[0] - 0.5) < 1e-9
    assert max(E) - min(E) < 0.01          # bounded, not drifting


def test_period_recovery():
    n = int(2 * math.pi / 0.001)
    qs, ps = velocity_verlet(lambda q: -q, 1.0, 0.0, 1.0, 0.001, n)
    assert abs(qs[-1][0] - 1.0) < 1e-3     # back to q(0) after one period


def test_leapfrog_equals_verlet_unit_mass():
    q1, p1 = velocity_verlet(lambda q: -q, 1.0, 0.0, 1.0, 0.01, 100)
    q2, v2 = leapfrog(lambda q: -q, 1.0, 0.0, 0.01, 100)
    assert max(abs(q1[i][0] - q2[i][0]) for i in range(101)) < 1e-12


def test_circular_orbit_radius_stable():
    def grav(q):
        x, y = q
        r = math.sqrt(x * x + y * y)
        return [-x / r ** 3, -y / r ** 3]
    qs, ps = velocity_verlet(grav, [1.0, 0.0], [0.0, 1.0], 1.0, 0.001, 10000)
    radii = [math.sqrt(q[0] ** 2 + q[1] ** 2) for q in qs]
    assert max(radii) - min(radii) < 0.01   # near-circular, radius preserved


def test_time_reversibility():
    # Integrate forward then backward (negate momentum) returns to start.
    qs, ps = velocity_verlet(lambda q: -q, 1.0, 0.3, 1.0, 0.01, 500)
    q_end, p_end = qs[-1][0], ps[-1][0]
    qb, pb = velocity_verlet(lambda q: -q, q_end, -p_end, 1.0, 0.01, 500)
    assert abs(qb[-1][0] - 1.0) < 1e-8
    assert abs(-pb[-1][0] - 0.3) < 1e-8


def test_validation():
    with pytest.raises(ValueError):
        velocity_verlet(lambda q: -q, 1.0, 0.0, 1.0, 0.1, 0)
