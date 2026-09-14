"""Total-variation denoising (Condat's fast 1-D algorithm).

Total-variation denoising recovers a piecewise-constant signal from a noisy one by
minimizing ``(1/2) sum (x - y)^2 + lam * sum |x_{k+1} - x_k|``. The L1 penalty on
successive differences flattens noise into constant runs while leaving sharp jumps
intact -- unlike a linear smoother, which rounds every edge. It is the standard
edge-preserving denoiser for step-like signals (piecewise-constant regimes, blocky
images row by row). Condat's algorithm solves the exact minimizer in a single ``O(n)``
forward pass. Pure standard library.
"""


def tv_denoise(y, lam):
    """Total-variation denoise ``y`` with regularization weight ``lam`` (Condat's method).

    Returns the exact minimizer of ``(1/2) sum (x - y)^2 + lam * sum |x_{k+1} - x_k|``.
    Larger ``lam`` produces flatter output with fewer, larger jumps; ``lam = 0`` returns
    ``y`` unchanged. The result is piecewise constant. ``O(n)`` time, ``O(1)`` extra
    state beyond the output.
    """
    n = len(y)
    if n == 0:
        raise ValueError("input must be non-empty")
    if lam < 0.0:
        raise ValueError("lam must be non-negative")
    x = [0.0] * n
    if lam == 0.0:
        return [float(v) for v in y]
    if n == 1:
        return [float(y[0])]

    k = 0            # current sample index being decided
    k0 = 0           # start of the current segment
    kplus = 0        # index of the last positive-jump candidate
    kminus = 0       # index of the last negative-jump candidate
    vmin = y[0] - lam
    vmax = y[0] + lam
    umin = lam
    umax = -lam

    while True:
        if k == n - 1:
            # Last sample: close out the current segment.
            if umin < 0.0:
                # negative jump: emit vmin up to kminus, restart
                for i in range(k0, kminus + 1):
                    x[i] = vmin
                k = k0 = kminus = kplus = kminus + 1
                vmin = y[k]
                vmax = y[k] + 2.0 * lam
                umin = lam
                umax = -lam
            elif umax > 0.0:
                # positive jump: emit vmax up to kplus, restart
                for i in range(k0, kplus + 1):
                    x[i] = vmax
                k = k0 = kplus = kminus = kplus + 1
                vmin = y[k] - 2.0 * lam
                vmax = y[k]
                umin = lam
                umax = -lam
            else:
                # no jump: level the tail
                vmin += umin / (k - k0 + 1)
                for i in range(k0, n):
                    x[i] = vmin
                break
            continue

        umin += y[k + 1] - vmin
        umax += y[k + 1] - vmax
        if umin < -lam:
            # negative jump confirmed at kminus
            for i in range(k0, kminus + 1):
                x[i] = vmin
            k = k0 = kplus = kminus = kminus + 1
            vmin = y[k]
            vmax = y[k] + 2.0 * lam
            umin = lam
            umax = -lam
        elif umax > lam:
            # positive jump confirmed at kplus
            for i in range(k0, kplus + 1):
                x[i] = vmax
            k = k0 = kplus = kminus = kplus + 1
            vmin = y[k] - 2.0 * lam
            vmax = y[k]
            umin = lam
            umax = -lam
        else:
            # continue the current segment; update running bounds
            k += 1
            if umin >= lam:
                kminus = k
                vmin += (umin - lam) / (kminus - k0 + 1)
                umin = lam
            if umax <= -lam:
                kplus = k
                vmax += (umax + lam) / (kplus - k0 + 1)
                umax = -lam
    return x


def tv_total_variation(x):
    """Total variation ``sum |x_{k+1} - x_k|`` of a sequence.

    The quantity the L1 penalty in :func:`tv_denoise` shrinks. Useful for confirming a
    denoised signal is flatter (lower total variation) than its noisy input.
    """
    return sum(abs(x[k + 1] - x[k]) for k in range(len(x) - 1))
