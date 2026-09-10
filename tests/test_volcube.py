"""Swaption vol cube: SABR-per-node with expiry/tenor variance interpolation."""

import pytest

from quantforge import VolCube
from quantforge.sabr import sabr_vol


TRUE = {
    (1, 5): (0.02, 0.5, -0.3, 0.4),
    (1, 10): (0.022, 0.5, -0.25, 0.45),
    (5, 5): (0.025, 0.5, -0.2, 0.35),
    (5, 10): (0.028, 0.5, -0.15, 0.3),
}
KS = [0.02, 0.03, 0.04, 0.05, 0.06]
F = 0.04


def _market():
    m = {}
    for (e, T), (al, be, rho, nu) in TRUE.items():
        vols = [sabr_vol(F, K, e, al, be, rho, nu) for K in KS]
        m[(e, T)] = (F, KS, vols)
    return m


def _cube():
    return VolCube.fit(_market(), beta=0.5)


@pytest.mark.parametrize("node", list(TRUE))
def test_node_reprices_its_smile(node):
    cube = _cube()
    e, T = node
    for K in KS:
        ref = sabr_vol(F, K, e, *TRUE[node])
        assert cube.vol(e, T, K) == pytest.approx(ref, abs=1e-3)


def test_interpolated_atm_vol_between_nodes():
    cube = _cube()
    v = cube.vol(3, 7, 0.04)
    v_lo = cube.vol(1, 5, 0.04)
    v_hi = cube.vol(5, 10, 0.04)
    assert min(v_lo, v_hi) - 0.02 < v < max(v_lo, v_hi) + 0.02


def test_expiry_interpolation_in_range():
    cube = _cube()
    v15 = cube.vol(1, 5, 0.04)
    v55 = cube.vol(5, 5, 0.04)
    mid = cube.vol(3, 5, 0.04)
    assert min(v15, v55) - 0.01 < mid < max(v15, v55) + 0.01


def test_smile_shape_preserved_at_node():
    cube = _cube()
    # Downward-then-up SABR smile: wings above the minimum.
    vols = [cube.vol(1, 5, K) for K in KS]
    assert vols[0] > min(vols) and vols[-1] > min(vols)


def test_node_smile_helper_matches_vol():
    cube = _cube()
    sm = cube.node_smile(1, 5, KS)
    for (K, iv), K2 in zip(sm, KS):
        assert K == K2
        assert iv == pytest.approx(cube.vol(1, 5, K), abs=1e-9)


def test_extrapolation_flat_beyond_grid():
    cube = _cube()
    # Beyond the tenor grid the nearest node is used (flat extrapolation).
    assert cube.vol(1, 20, 0.04) == pytest.approx(cube.vol(1, 10, 0.04), abs=1e-9)
