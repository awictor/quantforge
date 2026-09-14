"""Error-detecting and -correcting codes: Hamming(7,4) and the Luhn checksum.

The Hamming(7,4) code adds three parity bits to four data bits so that any single-bit
error is detected *and* corrected -- the textbook introduction to error-correcting codes.
The Luhn algorithm is the check-digit scheme guarding credit-card and many ID numbers,
catching all single-digit errors and most adjacent transpositions. Pure standard library.
"""


def hamming74_encode(bits):
    """Encode 4 data bits into a 7-bit Hamming(7,4) codeword.

    ``bits`` is a length-4 sequence of 0/1 (``d1 d2 d3 d4``). Returns a length-7 list
    ``[p1, p2, d1, p3, d2, d3, d4]`` where the parity bits ``p1, p2, p3`` cover the
    standard bit-position groups. Any single bit flip in the result is later correctable.
    """
    if len(bits) != 4 or any(b not in (0, 1) for b in bits):
        raise ValueError("bits must be four 0/1 values")
    d1, d2, d3, d4 = bits
    p1 = d1 ^ d2 ^ d4
    p2 = d1 ^ d3 ^ d4
    p3 = d2 ^ d3 ^ d4
    return [p1, p2, d1, p3, d2, d3, d4]


def hamming74_decode(code):
    """Decode a 7-bit Hamming(7,4) codeword, correcting any single-bit error.

    Computes the 3-bit syndrome; if non-zero it gives the 1-indexed position of the
    flipped bit, which is corrected before extracting the data. Returns
    ``(data_bits, error_position)`` where ``error_position`` is ``0`` if no error was
    found or the 1-indexed position that was corrected.
    """
    if len(code) != 7 or any(b not in (0, 1) for b in code):
        raise ValueError("code must be seven 0/1 values")
    c = list(code)
    # Positions are 1-indexed: parity bits at 1,2,4; data at 3,5,6,7.
    s1 = c[0] ^ c[2] ^ c[4] ^ c[6]      # covers positions 1,3,5,7
    s2 = c[1] ^ c[2] ^ c[5] ^ c[6]      # covers positions 2,3,6,7
    s3 = c[3] ^ c[4] ^ c[5] ^ c[6]      # covers positions 4,5,6,7
    syndrome = s1 * 1 + s2 * 2 + s3 * 4
    if syndrome != 0:
        c[syndrome - 1] ^= 1            # flip the erroneous bit back
    data = [c[2], c[4], c[5], c[6]]     # d1, d2, d3, d4
    return data, syndrome


def luhn_checksum(digits):
    """Luhn checksum of a digit sequence: ``0`` iff the number is valid.

    ``digits`` is a string or list of decimal digits *including* the trailing check digit.
    Doubles every second digit from the right, sums the digits of the results, and returns
    the total modulo 10. A valid Luhn number gives ``0``.
    """
    ds = [int(d) for d in digits]
    total = 0
    for i, d in enumerate(reversed(ds)):
        if i % 2 == 1:
            d *= 2
            if d > 9:
                d -= 9
        total += d
    return total % 10


def luhn_check_digit(digits):
    """Compute the Luhn check digit to append to a payload ``digits`` (without one).

    Returns the single digit ``0-9`` that makes ``digits + [check]`` a valid Luhn number
    (checksum zero).
    """
    ds = [int(d) for d in digits]
    total = 0
    # The appended check digit sits at position 0 from the right (not doubled), so the
    # payload's rightmost digit is doubled.
    for i, d in enumerate(reversed(ds)):
        if i % 2 == 0:
            d *= 2
            if d > 9:
                d -= 9
        total += d
    return (10 - total % 10) % 10
