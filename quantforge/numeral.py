"""Numeral systems: arbitrary-base conversion and Roman numerals.

Convert non-negative integers to and from any base ``2..36`` (digits ``0-9a-z``), and
between integers and Roman numerals. All conversions round-trip exactly. Pure standard
library.
"""

_DIGITS = "0123456789abcdefghijklmnopqrstuvwxyz"

_ROMAN = [
    (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
    (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
    (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I"),
]
_ROMAN_VALUES = {"I": 1, "V": 5, "X": 10, "L": 50, "C": 100, "D": 500, "M": 1000}


def to_base(n, base):
    """Represent a non-negative integer ``n`` in ``base`` (2..36) as a digit string.

    Digits above 9 use lowercase letters (``a`` = 10 .. ``z`` = 35). ``to_base(0, b)`` is
    ``"0"``.
    """
    if n < 0:
        raise ValueError("n must be non-negative")
    if not (2 <= base <= 36):
        raise ValueError("base must be in 2..36")
    if n == 0:
        return "0"
    out = []
    while n:
        out.append(_DIGITS[n % base])
        n //= base
    return "".join(reversed(out))


def from_base(s, base):
    """Parse a digit string ``s`` in ``base`` (2..36) back to a non-negative integer.

    Case-insensitive; raises on a digit that is out of range for the base.
    """
    if not (2 <= base <= 36):
        raise ValueError("base must be in 2..36")
    s = s.strip().lower()
    if not s:
        raise ValueError("empty string")
    value = 0
    for ch in s:
        d = _DIGITS.find(ch)
        if d < 0 or d >= base:
            raise ValueError("invalid digit %r for base %d" % (ch, base))
        value = value * base + d
    return value


def to_roman(n):
    """Convert an integer in ``1..3999`` to its Roman-numeral string (subtractive form)."""
    if not (1 <= n <= 3999):
        raise ValueError("Roman numerals cover 1..3999")
    out = []
    for value, sym in _ROMAN:
        while n >= value:
            out.append(sym)
            n -= value
    return "".join(out)


def from_roman(s):
    """Convert a Roman-numeral string to an integer (subtractive notation).

    Case-insensitive. Raises on characters that are not Roman digits.
    """
    s = s.strip().upper()
    if not s:
        raise ValueError("empty string")
    total = 0
    prev = 0
    for ch in reversed(s):
        if ch not in _ROMAN_VALUES:
            raise ValueError("invalid Roman digit %r" % (ch,))
        v = _ROMAN_VALUES[ch]
        if v < prev:
            total -= v
        else:
            total += v
            prev = v
    return total
