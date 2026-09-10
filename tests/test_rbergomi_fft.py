"""rBergomi FFT Volterra convolution matches the direct double loop."""

import pytest

from quantforge.rbergomi import rbergomi_paths


@pytest.mark.parametrize("n_steps", [32, 64, 128, 256])
def test_fft_matches_direct(n_steps):
    fast = rbergomi_paths(100, 1.0, 0.04, 1.5, 0.1, -0.7, 0.03,
                          n_steps=n_steps, n_paths=1500, seed=1, fast=True)
    direct = rbergomi_paths(100, 1.0, 0.04, 1.5, 0.1, -0.7, 0.03,
                            n_steps=n_steps, n_paths=1500, seed=1, fast=False)
    assert len(fast) == len(direct)
    assert max(abs(a - b) for a, b in zip(fast, direct)) < 1e-9


def test_auto_matches_forced_modes():
    # auto uses direct below 200 steps, FFT at/above -- both must equal an
    # explicit choice on the same seed.
    for n_steps, forced in [(100, False), (256, True)]:
        auto = rbergomi_paths(100, 1.0, 0.04, 1.5, 0.1, -0.7, 0.03,
                              n_steps=n_steps, n_paths=1000, seed=3)
        exp = rbergomi_paths(100, 1.0, 0.04, 1.5, 0.1, -0.7, 0.03,
                             n_steps=n_steps, n_paths=1000, seed=3, fast=forced)
        assert max(abs(a - b) for a, b in zip(auto, exp)) < 1e-9


def test_single_step_edge_case():
    # n_steps = 1 has no far cells; FFT path must degrade gracefully to direct.
    fast = rbergomi_paths(100, 0.5, 0.04, 1.0, 0.2, -0.5, 0.0,
                          n_steps=1, n_paths=500, seed=2, fast=True)
    direct = rbergomi_paths(100, 0.5, 0.04, 1.0, 0.2, -0.5, 0.0,
                            n_steps=1, n_paths=500, seed=2, fast=False)
    assert max(abs(a - b) for a, b in zip(fast, direct)) < 1e-12
