"""Forward-start variance-swap fair strike from two smiles."""

import pytest

from quantforge import (
    forward_variance_swap_from_smile as fwd,
    variance_swap_from_smile as vs,
)


S0, R = 100.0, 0.03


def _flat(sig):
    return lambda K: sig


def test_flat_term_structure_gives_variance():
    kf = fwd(S0, 0.5, 1.5, R, _flat(0.2), _flat(0.2))
    # Same finite-strip truncation as a spot var swap at sigma=0.2.
    assert kf == pytest.approx(vs(S0, 1.0, R, _flat(0.2)), abs=3e-3)


def test_t1_zero_is_spot_variance_swap():
    kf = fwd(S0, 0.0, 1.0, R, _flat(0.2), _flat(0.2))
    assert kf == pytest.approx(vs(S0, 1.0, R, _flat(0.2)), abs=1e-9)


def test_upward_term_structure_forward_above_front():
    kf = fwd(S0, 0.5, 1.5, R, _flat(0.18), _flat(0.22))
    front = vs(S0, 0.5, R, _flat(0.18))
    assert kf > front


def test_additivity_identity():
    # K_fwd (t2-t1) = K_var(t2) t2 - K_var(t1) t1.
    t1, t2 = 0.5, 2.0
    kf = fwd(S0, t1, t2, R, _flat(0.19), _flat(0.23))
    kv1 = vs(S0, t1, R, _flat(0.19))
    kv2 = vs(S0, t2, R, _flat(0.23))
    assert kf * (t2 - t1) == pytest.approx(kv2 * t2 - kv1 * t1, abs=1e-6)


def test_bad_times_raise():
    with pytest.raises(ValueError):
        fwd(S0, 1.5, 0.5, R, _flat(0.2), _flat(0.2))   # t1 > t2
