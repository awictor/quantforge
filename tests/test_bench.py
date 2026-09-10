"""Smoke test: the benchmark harness runs end-to-end on a tiny sample.

Keeps the benchmark importable and correct without timing assertions (which
would be machine-dependent and flaky in CI).
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "benchmarks"))

import bench  # noqa: E402


def test_sample_inputs_are_valid():
    rows = bench._sample_inputs(50, seed=1)
    assert len(rows) == 50
    for S, K, t, r, sigma in rows:
        assert S > 0 and K > 0 and t > 0 and sigma > 0


def test_bench_runs_small(capsys):
    bench.run(n=200, compare_vollib=False)
    out = capsys.readouterr().out
    assert "ops/sec" in out
    assert "round-trip" in out
    assert "solver iterations" in out


def test_newton_beats_bisection_iterations(capsys):
    # The Newton + Corrado-Miller solver should need far fewer iterations than
    # pure bisection; print both and check Newton is well under bisection.
    bench._iv_solver_iterations(bench._sample_inputs(2000, seed=3))
    out = capsys.readouterr().out
    lines = [ln for ln in out.splitlines() if "  " in ln and "." in ln]
    newton = float([ln for ln in lines if "Newton" in ln][0].split()[-1])
    bisect = float([ln for ln in lines if "bisection" in ln][0].split()[-1])
    assert newton < bisect / 2


def test_main_entrypoint(capsys):
    bench.main(["--n", "100"])
    assert "QuantForge benchmark" in capsys.readouterr().out
