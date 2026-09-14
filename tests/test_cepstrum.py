"""Real cepstrum, power cepstrum, quefrency (echo/pitch) detection."""

import random

import pytest

from quantforge import real_cepstrum, power_cepstrum, fundamental_quefrency


def test_length_is_next_power_of_two():
    x = [1.0] * 100          # 100 -> pads to 128
    assert len(real_cepstrum(x)) == 128
    x2 = [1.0] * 256
    assert len(real_cepstrum(x2)) == 256


def test_echo_delay_recovered():
    # y = x + 0.6 * x delayed by D  ->  cepstral peak at quefrency D
    N, D = 1024, 50
    rng = random.Random(42)
    base = [rng.gauss(0.0, 1.0) for _ in range(N)]
    y = [base[n] + (0.6 * base[n - D] if n >= D else 0.0) for n in range(N)]
    q, _ = fundamental_quefrency(y, min_quefrency=10, max_quefrency=300)
    assert q == D


def test_pitch_period_of_impulse_train():
    # voiced-speech model: a periodic impulse train -> peak at the period
    N, P = 1024, 64
    train = [1.0 if n % P == 0 else 0.0 for n in range(N)]
    x = [train[n] + 0.5 * train[n - 1] + 0.3 * train[n - 2] if n >= 2 else train[n]
         for n in range(N)]
    q, _ = fundamental_quefrency(x, min_quefrency=10, max_quefrency=300)
    assert q == P


def test_power_cepstrum_non_negative():
    N = 256
    rng = random.Random(1)
    x = [rng.gauss(0.0, 1.0) for _ in range(N)]
    assert min(power_cepstrum(x)) >= 0.0


def test_search_band_respected():
    N, D = 1024, 50
    rng = random.Random(7)
    base = [rng.gauss(0.0, 1.0) for _ in range(N)]
    y = [base[n] + 0.6 * base[n - D] if n >= D else base[n] for n in range(N)]
    # exclude the true peak; the reported quefrency must fall inside the band
    q, _ = fundamental_quefrency(y, min_quefrency=60, max_quefrency=120)
    assert 60 <= q <= 120


def test_validation():
    with pytest.raises(ValueError):
        real_cepstrum([])
    with pytest.raises(ValueError):
        power_cepstrum([])
    with pytest.raises(ValueError):
        fundamental_quefrency([1.0] * 64, min_quefrency=40, max_quefrency=20)
