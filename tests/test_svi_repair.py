"""Tests for single-slice SVI butterfly-arbitrage repair."""

import pytest

from quantforge import (
    SVIParams, svi_repair_butterfly, svi_is_butterfly_free,
    svi_butterfly_arbitrage,
)


def test_repairs_pathological_slice():
    bad = SVIParams(a=0.001, b=0.9, rho=-0.9, m=0.0, s=0.02)
    assert not svi_is_butterfly_free(bad)
    fixed = svi_repair_butterfly(bad)
    assert svi_is_butterfly_free(fixed)


def test_repair_shrinks_wings():
    bad = SVIParams(a=0.001, b=0.9, rho=-0.9, m=0.0, s=0.02)
    fixed = svi_repair_butterfly(bad)
    assert fixed.b < bad.b


def test_repair_preserves_other_params():
    bad = SVIParams(a=0.001, b=0.9, rho=-0.9, m=0.05, s=0.02)
    fixed = svi_repair_butterfly(bad)
    assert fixed.a == bad.a
    assert fixed.rho == bad.rho
    assert fixed.m == bad.m
    assert fixed.s == bad.s


def test_clean_slice_returned_unchanged():
    good = SVIParams(a=0.04, b=0.2, rho=-0.3, m=0.0, s=0.3)
    assert svi_repair_butterfly(good) is good


def test_repaired_slice_has_no_violations():
    bad = SVIParams(a=0.002, b=0.8, rho=-0.7, m=0.0, s=0.05)
    fixed = svi_repair_butterfly(bad)
    assert svi_butterfly_arbitrage(fixed) == []
