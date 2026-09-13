"""Multiple-hypothesis-testing corrections."""

import pytest

from quantforge import bonferroni, holm, benjamini_hochberg, benjamini_yekutieli

ALL = [bonferroni, holm, benjamini_hochberg, benjamini_yekutieli]
P = [0.01, 0.02, 0.03, 0.04, 0.05]


def test_bonferroni_values():
    assert bonferroni(P) == [0.05, 0.10, 0.15, 0.20, 0.25]


def test_bh_known_example():
    # statsmodels fdr_bh on this input gives all 0.05.
    adj = benjamini_hochberg(P)
    assert all(abs(a - 0.05) < 1e-12 for a in adj)


def test_by_is_bh_times_harmonic():
    h = sum(1.0 / i for i in range(1, len(P) + 1))
    bh = benjamini_hochberg(P)
    by = benjamini_yekutieli(P)
    for a, b in zip(by, bh):
        assert abs(a - min(1.0, b * h)) < 1e-12


def test_adjusted_at_least_raw():
    for f in ALL:
        adj = f(P)
        assert all(adj[i] >= P[i] - 1e-12 for i in range(len(P)))
        assert all(0.0 <= a <= 1.0 for a in adj)


def test_conservativeness_ordering():
    b = bonferroni(P)
    h = holm(P)
    bh = benjamini_hochberg(P)
    by = benjamini_yekutieli(P)
    # Holm dominates Bonferroni; BY is more conservative than BH.
    assert all(b[i] >= h[i] - 1e-12 for i in range(len(P)))
    assert all(by[i] >= bh[i] - 1e-12 for i in range(len(P)))
    # Family-wise (Holm) at least as large as FDR (BH).
    assert all(h[i] >= bh[i] - 1e-12 for i in range(len(P)))


def test_order_preserved_for_unsorted_input():
    pu = [0.04, 0.01, 0.05, 0.02, 0.03]
    adj = benjamini_hochberg(pu)
    # Reordering the input reorders the output the same way.
    ref = benjamini_hochberg(sorted(pu))
    idx = sorted(range(len(pu)), key=lambda i: pu[i])
    for rank, i in enumerate(idx):
        assert abs(adj[i] - ref[rank]) < 1e-12


def test_monotone_in_rank():
    # Sorted p-values give non-decreasing adjusted values for step methods.
    for f in (holm, benjamini_hochberg, benjamini_yekutieli, bonferroni):
        adj = f(P)                       # P already ascending
        assert all(adj[i] <= adj[i + 1] + 1e-12 for i in range(len(P) - 1))


def test_single_pvalue_unchanged():
    for f in ALL:
        assert abs(f([0.037])[0] - 0.037) < 1e-12


def test_all_tiny_all_reject():
    tiny = [1e-6, 2e-6, 3e-6]
    for f in ALL:
        assert all(a < 0.05 for a in f(tiny))


def test_validation():
    for f in ALL:
        with pytest.raises(ValueError):
            f([])
        with pytest.raises(ValueError):
            f([0.5, 1.5])
