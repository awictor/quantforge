"""Throughput and accuracy benchmarks for QuantForge.

Pure standard library (``time.perf_counter``) — no pytest, no third-party
timing harness. Run:

    python benchmarks/bench.py               # QuantForge only
    python benchmarks/bench.py --vollib       # also compare to py_vollib if installed
    python benchmarks/bench.py --n 200000     # set the sample size

Prints a table of operations/second and (when py_vollib is present) the max
absolute pricing difference, so the README's speed/accuracy claims are
reproducible on any machine.
"""

import argparse
import math
import random
import statistics
import time

from quantforge import (
    call_price, greeks, implied_volatility, american_price, OptionType,
)


def _sample_inputs(n, seed=0):
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        S = rng.uniform(50, 150)
        K = rng.uniform(50, 150)
        t = rng.uniform(0.05, 2.0)
        r = rng.uniform(0.0, 0.08)
        sigma = rng.uniform(0.05, 0.8)
        out.append((S, K, t, r, sigma))
    return out


def _timed(label, fn, inputs):
    start = time.perf_counter()
    for args in inputs:
        fn(*args)
    elapsed = time.perf_counter() - start
    ops = len(inputs) / elapsed if elapsed > 0 else float("inf")
    print(f"  {label:<28}{ops:>14,.0f} ops/sec  ({elapsed*1e6/len(inputs):.2f} us/op)")
    return ops


def run(n, compare_vollib):
    inputs = _sample_inputs(n)
    print(f"QuantForge benchmark  (n={n:,} contracts)\n")

    print("Throughput:")
    _timed("call price (BSM)", lambda S, K, t, r, s: call_price(S, K, t, r, s), inputs)
    _timed("full greeks (6 outputs)",
           lambda S, K, t, r, s: greeks(S, K, t, r, s, OptionType.CALL), inputs)

    # Implied vol on prices we generated, so every solve is well-posed.
    iv_inputs = []
    for (S, K, t, r, sigma) in inputs:
        p = call_price(S, K, t, r, sigma)
        iv_inputs.append((p, S, K, t, r))
    _timed("implied vol (Newton+bisect)",
           lambda p, S, K, t, r: implied_volatility(p, S, K, t, r, OptionType.CALL),
           iv_inputs)

    # American pricing is O(steps^2); use a smaller sample.
    am_inputs = inputs[: max(1, n // 100)]
    _timed("american (200-step tree)",
           lambda S, K, t, r, s: american_price(S, K, t, r, s, OptionType.CALL, steps=200),
           am_inputs)

    # Accuracy: implied vol round-trip error.
    errs = []
    for (S, K, t, r, sigma) in inputs[: min(n, 20000)]:
        p = call_price(S, K, t, r, sigma)
        try:
            iv = implied_volatility(p, S, K, t, r, OptionType.CALL)
            errs.append(abs(call_price(S, K, t, r, iv) - p))
        except ValueError:
            pass
    print("\nAccuracy (implied-vol round-trip price error):")
    print(f"  max  {max(errs):.2e}")
    print(f"  mean {statistics.fmean(errs):.2e}")

    if compare_vollib:
        _compare_vollib(inputs)


def _compare_vollib(inputs):
    try:
        from py_vollib.black_scholes import black_scholes as bs_vollib
    except ImportError:
        print("\npy_vollib not installed; skipping accuracy comparison.")
        print("  pip install py_vollib   # to enable")
        return
    max_diff = 0.0
    for (S, K, t, r, sigma) in inputs[: min(len(inputs), 20000)]:
        ours = call_price(S, K, t, r, sigma)
        theirs = bs_vollib("c", S, K, t, r, sigma)
        max_diff = max(max_diff, abs(ours - theirs))
    print(f"\nvs py_vollib (call price):")
    print(f"  max abs diff  {max_diff:.2e}")


def main(argv=None):
    ap = argparse.ArgumentParser(description="QuantForge benchmarks")
    ap.add_argument("--n", type=int, default=100_000, help="number of contracts")
    ap.add_argument("--vollib", action="store_true", help="compare to py_vollib")
    args = ap.parse_args(argv)
    run(args.n, args.vollib)


if __name__ == "__main__":
    main()
