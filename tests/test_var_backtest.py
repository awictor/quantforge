"""VaR / ES backtests: Kupiec, Christoffersen, Acerbi-Szekely."""

import math
import random

import pytest

from quantforge import (kupiec_pof, christoffersen_independence,
                        christoffersen_cc, acerbi_szekely_es)
from quantforge.mathfns import norm_ppf

ALPHA = 0.01
Z = norm_ppf(1 - ALPHA)
ES = math.exp(-0.5 * Z * Z) / math.sqrt(2 * math.pi) / ALPHA


def _normal_losses(n, seed):
    rng = random.Random(seed)
    return [rng.gauss(0, 1) for _ in range(n)]


def test_kupiec_accepts_correct_model():
    losses = _normal_losses(5000, 1)
    _, p = kupiec_pof(losses, [Z] * 5000, ALPHA)
    assert p > 0.05


def test_kupiec_rejects_too_many_exceptions():
    losses = _normal_losses(5000, 1)
    _, p = kupiec_pof(losses, [1.5] * 5000, ALPHA)   # VaR far too low
    assert p < 0.01


def test_christoffersen_independence_accepts_iid():
    losses = _normal_losses(5000, 2)
    _, p = christoffersen_independence(losses, [Z] * 5000)
    assert p > 0.05


def test_christoffersen_cc_accepts_correct():
    losses = _normal_losses(5000, 3)
    _, p = christoffersen_cc(losses, [Z] * 5000, ALPHA)
    assert p > 0.05


def test_acerbi_szekely_near_zero_when_calibrated():
    losses = _normal_losses(5000, 1)
    z = acerbi_szekely_es(losses, [Z] * 5000, [ES] * 5000, ALPHA)
    assert abs(z) < 0.2


def test_acerbi_szekely_positive_when_es_understated():
    losses = _normal_losses(5000, 1)
    z = acerbi_szekely_es(losses, [Z] * 5000, [ES * 0.5] * 5000, ALPHA)
    assert z > 0.3


def test_acerbi_szekely_negative_when_es_overstated():
    losses = _normal_losses(5000, 1)
    z = acerbi_szekely_es(losses, [Z] * 5000, [ES * 2.0] * 5000, ALPHA)
    assert z < 0.0


def test_validation():
    with pytest.raises(ValueError):
        kupiec_pof([1.0], [2.0], alpha=1.5)
    with pytest.raises(ValueError):
        kupiec_pof([], [], ALPHA)
    with pytest.raises(ValueError):
        acerbi_szekely_es([1.0], [2.0], [0.0], ALPHA)   # non-positive ES
