"""Bakshi-Kapadia-Madan risk-neutral moments from a smile."""

import math

import pytest

from quantforge import bkm_moments_from_smile, skew_swap_from_smile


S0, T, R = 100.0, 0.5, 0.03


def _flat(K):
    return 0.2


def test_flat_smile_is_symmetric():
    var, skew, exk = bkm_moments_from_smile(S0, T, R, _flat,
                                            n_strikes=401, width=8.0)
    assert var == pytest.approx(0.04 * T, abs=1e-3)   # sigma^2 t
    assert skew == pytest.approx(0.0, abs=1e-2)
    assert exk == pytest.approx(0.0, abs=5e-2)


def test_downward_skew_negative_skewness_fat_tails():
    def down(K):
        return max(0.05, 0.2 + 0.2 * math.log(S0 / K))
    var, skew, exk = bkm_moments_from_smile(S0, T, R, down,
                                            n_strikes=401, width=8.0)
    assert skew < 0.0
    assert exk > 0.0


def test_upward_skew_positive_skewness():
    def up(K):
        return max(0.05, 0.2 - 0.2 * math.log(S0 / K))
    _var, skew, _exk = bkm_moments_from_smile(S0, T, R, up,
                                              n_strikes=401, width=8.0)
    assert skew > 0.0


def test_skew_swap_matches_moment_skew():
    def down(K):
        return max(0.05, 0.2 + 0.15 * math.log(S0 / K))
    ss = skew_swap_from_smile(S0, T, R, down, n_strikes=401, width=8.0)
    _var, skew, _exk = bkm_moments_from_smile(S0, T, R, down,
                                              n_strikes=401, width=8.0)
    assert ss == pytest.approx(skew, abs=1e-9)


def test_steeper_skew_more_negative():
    def sk(slope):
        return lambda K: max(0.05, 0.2 + slope * math.log(S0 / K))
    mild = skew_swap_from_smile(S0, T, R, sk(0.1), n_strikes=401, width=8.0)
    steep = skew_swap_from_smile(S0, T, R, sk(0.25), n_strikes=401, width=8.0)
    assert steep < mild


def test_bad_tenor_raises():
    with pytest.raises(ValueError):
        bkm_moments_from_smile(S0, 0.0, R, _flat)
