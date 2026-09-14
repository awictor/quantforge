"""Classic sequence algorithms: LIS, maximum subarray, longest run."""

import random

import pytest

from quantforge import (
    longest_increasing_subsequence,
    maximum_subarray,
    longest_run,
)


def _brute_lis_len(x, strict):
    n = len(x)
    if n == 0:
        return 0
    dp = [1] * n
    for i in range(n):
        for j in range(i):
            if (x[j] < x[i]) if strict else (x[j] <= x[i]):
                dp[i] = max(dp[i], dp[j] + 1)
    return max(dp)


def _is_subseq(sub, full):
    it = iter(full)
    return all(c in it for c in sub)


def test_lis_vs_brute():
    rng = random.Random(1)
    for _ in range(1000):
        x = [rng.randint(0, 15) for _ in range(rng.randint(0, 20))]
        for strict in (True, False):
            seq = longest_increasing_subsequence(x, strict=strict)
            assert len(seq) == _brute_lis_len(x, strict)
            assert _is_subseq(seq, x)
            for i in range(1, len(seq)):
                assert (seq[i] > seq[i - 1]) if strict else (seq[i] >= seq[i - 1])
    assert longest_increasing_subsequence([10, 9, 2, 5, 3, 7, 101, 18]) == [2, 3, 7, 18]


def test_kadane_vs_brute():
    def brute(x):
        best = x[0]
        for i in range(len(x)):
            s = 0
            for j in range(i, len(x)):
                s += x[j]
                best = max(best, s)
        return best

    rng = random.Random(2)
    for _ in range(1000):
        x = [rng.randint(-10, 10) for _ in range(rng.randint(1, 25))]
        s, lo, hi = maximum_subarray(x)
        assert s == brute(x)
        assert sum(x[lo:hi + 1]) == s
    assert maximum_subarray([-2, 1, -3, 4, -1, 2, 1, -5, 4])[0] == 6
    assert maximum_subarray([-5, -2, -8])[0] == -2


def test_longest_run():
    assert longest_run([1, 1, 2, 3, 3, 3, 1]) == (3, 3, 3)
    rng = random.Random(3)
    for _ in range(500):
        x = [rng.randint(0, 3) for _ in range(rng.randint(1, 20))]
        val, length, start = longest_run(x)
        assert all(x[start + i] == val for i in range(length))
        assert start == 0 or x[start - 1] != val
        assert start + length == len(x) or x[start + length] != val


def test_validation():
    with pytest.raises(ValueError):
        maximum_subarray([])
    with pytest.raises(ValueError):
        longest_run([])
