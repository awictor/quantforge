"""Sample-rate conversion: sinc interpolation, up/downsampling, rational resample."""

import math

import pytest

from quantforge import sinc_interp, upsample, downsample, resample_rational

NUMTAPS = 97
DELAY = (NUMTAPS - 1) // 2


def test_sinc_interp_exact_at_integers():
    x = [1.0, 3.0, -2.0, 5.0, 0.5, -1.0, 2.0, 4.0]
    back = sinc_interp(x, list(range(len(x))))
    assert max(abs(back[i] - x[i]) for i in range(len(x))) < 1e-9


def test_sinc_interp_reconstructs_tone():
    N, f = 64, 3
    tone = [math.cos(2 * math.pi * f * n / N) for n in range(N)]
    mid = sinc_interp(tone, [n + 0.5 for n in range(15, N - 15)])
    true = [math.cos(2 * math.pi * f * (n + 0.5) / N) for n in range(15, N - 15)]
    assert max(abs(mid[i] - true[i]) for i in range(len(mid))) < 0.02


def test_upsample_preserves_low_tone():
    N = 64
    lo = [math.cos(2 * math.pi * 2 * n / N) for n in range(N)]
    up = upsample(lo, 4, numtaps=NUMTAPS)
    assert len(up) == 4 * N
    # up[m] aligns to input continuous time (m - DELAY)/4
    errs = [abs(up[m] - math.cos(2 * math.pi * 2 * ((m - DELAY) / 4) / N))
            for m in range(120, len(up) - 40)]
    assert max(errs) < 0.01


def test_downsample_preserves_low_tone():
    N = 256
    lo = [math.cos(2 * math.pi * 8 * n / N) for n in range(N)]
    d = downsample(lo, 2, numtaps=NUMTAPS)
    assert len(d) == (N + 1) // 2
    errs = [abs(d[m] - math.cos(2 * math.pi * 8 * (2 * m - DELAY) / N))
            for m in range(40, len(d) - 5)]
    assert max(errs) < 0.01


def test_downsample_anti_aliases():
    # a tone above the new Nyquist (0.32 > 0.25 after /2) must be suppressed, not aliased
    N = 256
    hi = [math.cos(2 * math.pi * 0.32 * n) for n in range(N)]
    d = downsample(hi, 2, numtaps=NUMTAPS)
    rms = (sum(v * v for v in d[40:]) / len(d[40:])) ** 0.5
    assert rms < 0.05


def test_rational_resample_length_and_shape():
    N = 256
    lo = [math.cos(2 * math.pi * 8 * n / N) for n in range(N)]
    r = resample_rational(lo, 3, 2, numtaps=NUMTAPS)
    assert abs(len(r) - N * 3 // 2) <= 2
    errs = [abs(r[m] - math.cos(2 * math.pi * 8 * ((2 * m - DELAY) / 3) / N))
            for m in range(80, len(r) - 40)]
    assert max(errs) < 0.01


def test_factor_one_is_identity():
    x = [1.0, 2.0, 3.0, 4.0]
    assert upsample(x, 1) == x
    assert downsample(x, 1) == x
    assert resample_rational(x, 1, 1) == x


def test_validation():
    with pytest.raises(ValueError):
        upsample([], 2)
    with pytest.raises(ValueError):
        downsample([1.0, 2.0, 3.0], 0)
    with pytest.raises(ValueError):
        sinc_interp([], [0.0])
    with pytest.raises(ValueError):
        resample_rational([1.0, 2.0], 0, 2)
