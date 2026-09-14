"""Short-time Fourier transform, spectrogram, inverse STFT."""

import math

import pytest

from quantforge import stft, spectrogram, istft


def test_roundtrip_reconstruction():
    N = 2048
    x = [math.sin(2 * math.pi * 5 * n / 256) +
         0.3 * math.sin(2 * math.pi * 20 * n / 256) for n in range(N)]
    frames = stft(x, 256, hop=64, window="hann")
    xr = istft(frames, 256, hop=64, window="hann", length=N)
    # interior avoids partial window coverage at the very edges
    assert max(abs(xr[n] - x[n]) for n in range(256, N - 256)) < 1e-10


def test_tone_concentrates_at_its_bin():
    fs, k = 256, 10
    tone = [math.cos(2 * math.pi * k * n / fs) for n in range(1024)]
    S = spectrogram(tone, 256, hop=128, window="hann")
    mid = S[len(S) // 2]
    assert max(range(len(mid)), key=lambda b: mid[b]) == k


def test_chirp_peak_rises_over_time():
    T = 4096
    chirp = [math.sin(2 * math.pi * (2 + 40 * n / T) * n / T) for n in range(T)]
    S = spectrogram(chirp, 256, hop=128, window="hann")
    peaks = [max(range(len(col)), key=lambda b: col[b]) for col in S]
    early = sum(peaks[:5]) / 5
    late = sum(peaks[-5:]) / 5
    assert late > early


def test_spectrogram_dimensions():
    x = [math.sin(n / 3.0) for n in range(1024)]
    S = spectrogram(x, 256, hop=128, window="hann")
    assert len(S[0]) == 256 // 2 + 1          # DC..Nyquist
    assert len(S) == len(stft(x, 256, hop=128))


def test_hop_controls_frame_count():
    x = [1.0] * 1024
    assert len(stft(x, 256, hop=256)) == 4     # non-overlapping
    assert len(stft(x, 256, hop=128)) == 8     # 50% overlap


def test_validation():
    with pytest.raises(ValueError):
        stft([1.0] * 100, 255)                  # frame_size not power of two
    with pytest.raises(ValueError):
        stft([], 256)                           # empty
    with pytest.raises(ValueError):
        stft([1.0] * 100, 256, hop=0)           # bad hop
    with pytest.raises(ValueError):
        istft([], 256)                          # no frames
