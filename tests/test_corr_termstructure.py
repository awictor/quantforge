"""Tests for the implied-correlation term structure."""

import pytest

from quantforge import (
    correlation_term_structure, index_vol_from_correlation,
)


W = [0.4, 0.35, 0.25]


def _member_curves():
    # Three members, each with a rising vol term structure across 3 expiries.
    return [
        [0.22, 0.24, 0.26],   # member 0
        [0.28, 0.30, 0.31],   # member 1
        [0.18, 0.19, 0.20],   # member 2
    ]


def test_recovers_constant_correlation():
    curves = _member_curves()
    expiries = [0.5, 1.0, 2.0]
    rho_true = 0.35
    # Build the index vol at each expiry from a constant correlation.
    index_curve = []
    for j in range(len(expiries)):
        vols_j = [curves[i][j] for i in range(len(W))]
        index_curve.append(index_vol_from_correlation(W, vols_j, rho_true))

    ts = correlation_term_structure(W, curves, index_curve, expiries)
    assert [e for e, _ in ts] == expiries
    for _, rho in ts:
        assert rho == pytest.approx(rho_true, abs=1e-9)


def test_recovers_varying_correlation():
    curves = _member_curves()
    expiries = [0.5, 1.0, 2.0]
    rhos_true = [0.2, 0.4, 0.6]   # rising correlation with maturity
    index_curve = []
    for j in range(len(expiries)):
        vols_j = [curves[i][j] for i in range(len(W))]
        index_curve.append(index_vol_from_correlation(W, vols_j, rhos_true[j]))

    ts = correlation_term_structure(W, curves, index_curve, expiries)
    for (_, rho), rt in zip(ts, rhos_true):
        assert rho == pytest.approx(rt, abs=1e-9)


def test_length_matches_expiries():
    curves = _member_curves()
    expiries = [0.5, 1.0, 2.0]
    index_curve = [0.2, 0.22, 0.24]
    ts = correlation_term_structure(W, curves, index_curve, expiries)
    assert len(ts) == 3


def test_rejects_mismatched_lengths():
    curves = _member_curves()
    with pytest.raises(ValueError):
        correlation_term_structure(W, curves, [0.2, 0.22], [0.5, 1.0, 2.0])
    with pytest.raises(ValueError):
        bad = [[0.2, 0.2], [0.3, 0.3, 0.3], [0.2, 0.2, 0.2]]
        correlation_term_structure(W, bad, [0.2, 0.22, 0.24], [0.5, 1.0, 2.0])
