"""Color-space conversions: RGB, HSV, HSL, and hex.

Conversions between the common color representations. RGB channels are floats in
``[0, 1]``; HSV/HSL hue is in degrees ``[0, 360)`` with saturation/value/lightness in
``[0, 1]``. HSV (hue-saturation-value) and HSL (hue-saturation-lightness) are the
cylindrical models used by color pickers; hex is the ``"#rrggbb"`` web form. All
conversions round-trip. Pure standard library.
"""


def rgb_to_hsv(r, g, b):
    """Convert RGB (each in ``[0, 1]``) to ``(hue_deg, saturation, value)``.

    Hue in ``[0, 360)`` (0 for gray), saturation and value in ``[0, 1]``.
    """
    mx = max(r, g, b)
    mn = min(r, g, b)
    delta = mx - mn
    if delta == 0:
        h = 0.0
    elif mx == r:
        h = 60.0 * (((g - b) / delta) % 6.0)
    elif mx == g:
        h = 60.0 * ((b - r) / delta + 2.0)
    else:
        h = 60.0 * ((r - g) / delta + 4.0)
    s = 0.0 if mx == 0 else delta / mx
    return h % 360.0, s, mx


def hsv_to_rgb(h, s, v):
    """Convert ``(hue_deg, saturation, value)`` to RGB (each in ``[0, 1]``)."""
    h = h % 360.0
    c = v * s
    x = c * (1.0 - abs((h / 60.0) % 2.0 - 1.0))
    m = v - c
    if h < 60:
        rp, gp, bp = c, x, 0.0
    elif h < 120:
        rp, gp, bp = x, c, 0.0
    elif h < 180:
        rp, gp, bp = 0.0, c, x
    elif h < 240:
        rp, gp, bp = 0.0, x, c
    elif h < 300:
        rp, gp, bp = x, 0.0, c
    else:
        rp, gp, bp = c, 0.0, x
    return rp + m, gp + m, bp + m


def rgb_to_hsl(r, g, b):
    """Convert RGB (each in ``[0, 1]``) to ``(hue_deg, saturation, lightness)``."""
    mx = max(r, g, b)
    mn = min(r, g, b)
    delta = mx - mn
    l = (mx + mn) / 2.0
    if delta == 0:
        return 0.0, 0.0, l
    s = delta / (1.0 - abs(2.0 * l - 1.0))
    if mx == r:
        h = 60.0 * (((g - b) / delta) % 6.0)
    elif mx == g:
        h = 60.0 * ((b - r) / delta + 2.0)
    else:
        h = 60.0 * ((r - g) / delta + 4.0)
    return h % 360.0, s, l


def hsl_to_rgb(h, s, l):
    """Convert ``(hue_deg, saturation, lightness)`` to RGB (each in ``[0, 1]``)."""
    h = h % 360.0
    c = (1.0 - abs(2.0 * l - 1.0)) * s
    x = c * (1.0 - abs((h / 60.0) % 2.0 - 1.0))
    m = l - c / 2.0
    if h < 60:
        rp, gp, bp = c, x, 0.0
    elif h < 120:
        rp, gp, bp = x, c, 0.0
    elif h < 180:
        rp, gp, bp = 0.0, c, x
    elif h < 240:
        rp, gp, bp = 0.0, x, c
    elif h < 300:
        rp, gp, bp = x, 0.0, c
    else:
        rp, gp, bp = c, 0.0, x
    return rp + m, gp + m, bp + m


def rgb_to_hex(r, g, b):
    """Convert RGB (each in ``[0, 1]``) to a ``"#rrggbb"`` hex string (rounded, clamped)."""
    def ch(v):
        return max(0, min(255, int(round(v * 255.0))))
    return "#%02x%02x%02x" % (ch(r), ch(g), ch(b))


def hex_to_rgb(code):
    """Convert a ``"#rrggbb"`` (or ``"rrggbb"``) hex string to RGB floats in ``[0, 1]``."""
    s = code.lstrip("#")
    if len(s) != 6:
        raise ValueError("hex code must be 6 digits (optionally prefixed with #)")
    try:
        r = int(s[0:2], 16)
        g = int(s[2:4], 16)
        b = int(s[4:6], 16)
    except ValueError:
        raise ValueError("invalid hex digits")
    return r / 255.0, g / 255.0, b / 255.0
