"""Latin hypercube sampling and uniformity diagnostics."""

import statistics

import pytest

from quantforge import latin_hypercube, maximin_lhs, l2_star_discrepancy
from quantforge.lhs import _lcg, _min_pairwise_distance


def test_stratification_one_point_per_bin():
    n, dim = 20, 3
    pts = latin_hypercube(n, dim, seed=42)
    assert len(pts) == n and all(len(p) == dim for p in pts)
    for d in range(dim):
        bins = sorted(int(p[d] * n) for p in pts)
        assert bins == list(range(n))


def test_centered_points_at_bin_centers():
    c = latin_hypercube(5, 2, seed=1, centered=True)
    for p in c:
        for x in p:
            # each coordinate is (integer + 0.5)/5
            assert abs((x * 5) - (round(x * 5 - 0.5) + 0.5)) < 1e-9


def test_lhs_more_uniform_than_random():
    lhs_disc = [l2_star_discrepancy(latin_hypercube(30, 2, seed=s)) for s in range(20)]

    def random_set(n, dim, seed):
        r = _lcg(seed)
        return [[r() for _ in range(dim)] for _ in range(n)]

    rnd_disc = [l2_star_discrepancy(random_set(30, 2, seed=s + 999)) for s in range(20)]
    assert statistics.mean(lhs_disc) < statistics.mean(rnd_disc)


def test_maximin_improves_spread():
    single = _min_pairwise_distance(latin_hypercube(15, 2, seed=7))
    mm = _min_pairwise_distance(maximin_lhs(15, 2, seed=7, tries=30))
    assert mm >= single


def test_discrepancy_positive():
    assert l2_star_discrepancy([[0.5]]) > 0.0
    assert l2_star_discrepancy(latin_hypercube(10, 3, seed=3)) > 0.0


def test_validation():
    with pytest.raises(ValueError):
        latin_hypercube(0, 2)
    with pytest.raises(ValueError):
        latin_hypercube(5, 0)
    with pytest.raises(ValueError):
        maximin_lhs(5, 2, tries=0)
    with pytest.raises(ValueError):
        l2_star_discrepancy([])
