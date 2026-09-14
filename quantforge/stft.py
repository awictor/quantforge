"""Short-time Fourier transform and spectrogram (time-frequency analysis).

A single FFT tells you *what* frequencies are in a signal but not *when* they occur. The
short-time Fourier transform slides a window along the signal and takes an FFT of each
overlapping frame, giving a two-dimensional picture -- frequency against time. Its
squared magnitude is the *spectrogram*, the standard view of speech, music, and any
signal whose spectrum changes over time. The transform is invertible by weighted
overlap-add, so you can filter in the time-frequency plane and reconstruct. Pure
standard library on top of the radix-2 FFT and the window functions.
"""

from .fft import fft, ifft
from .windows import apply_window, hann


def _is_pow2(n):
    return n >= 1 and (n & (n - 1)) == 0


def stft(x, frame_size, hop=None, window="hann"):
    """Short-time Fourier transform: a list of per-frame complex FFT spectra.

    Slides a ``frame_size`` window along ``x`` in steps of ``hop`` (default
    ``frame_size // 2``, i.e. 50% overlap), applies the ``window`` taper, and FFTs each
    frame. Returns a list of frames, each a length-``frame_size`` complex spectrum.
    ``frame_size`` must be a power of two (FFT constraint). Frames that run past the end
    are zero-padded.
    """
    n = len(x)
    if not _is_pow2(frame_size):
        raise ValueError("frame_size must be a power of two")
    if n == 0:
        raise ValueError("input must be non-empty")
    if hop is None:
        hop = frame_size // 2
    if hop < 1:
        raise ValueError("hop must be >= 1")
    frames = []
    start = 0
    while start < n:
        seg = list(x[start:start + frame_size])
        if len(seg) < frame_size:
            seg = seg + [0.0] * (frame_size - len(seg))
        frames.append(fft(apply_window(seg, window)))
        start += hop
    return frames


def spectrogram(x, frame_size, hop=None, window="hann"):
    """Spectrogram: per-frame power ``|STFT|^2`` over the lower (non-redundant) half.

    Returns a list of frames, each a list of ``frame_size // 2 + 1`` power values (DC up
    to Nyquist -- the rest of a real signal's spectrum is a mirror image). Each frame is
    a column in time; each entry a frequency bin.
    """
    half = frame_size // 2 + 1
    return [[abs(v) ** 2 for v in frame[:half]] for frame in stft(x, frame_size, hop, window)]


def istft(frames, frame_size, hop=None, window="hann", length=None):
    """Invert an :func:`stft` back to the time domain by weighted overlap-add.

    Applies the synthesis window to each inverse-FFT frame and overlap-adds, dividing by
    the overlap-added squared window so the reconstruction is exact wherever the window
    coverage is non-zero (independent of the constant-overlap-add condition). ``length``
    truncates the output; by default it keeps the full overlap-added span.
    """
    if not frames:
        raise ValueError("no frames")
    if not _is_pow2(frame_size):
        raise ValueError("frame_size must be a power of two")
    if hop is None:
        hop = frame_size // 2
    win = window if isinstance(window, list) else _named_window(window, frame_size)
    total = hop * (len(frames) - 1) + frame_size
    out = [0.0] * total
    norm = [0.0] * total
    for i, frame in enumerate(frames):
        time_frame = [v.real for v in ifft(frame)]
        start = i * hop
        for k in range(frame_size):
            out[start + k] += time_frame[k] * win[k]
            norm[start + k] += win[k] * win[k]
    result = [out[n] / norm[n] if norm[n] > 1e-12 else 0.0 for n in range(total)]
    if length is not None:
        result = result[:length]
    return result


def _named_window(name, n):
    from .windows import hamming, blackman, bartlett, rectangular
    table = {
        "hann": hann,
        "hamming": hamming,
        "blackman": blackman,
        "bartlett": bartlett,
        "rectangular": rectangular,
        "boxcar": rectangular,
    }
    if name not in table:
        raise ValueError("unknown window %r" % (name,))
    return table[name](n)
