"""Checksums and non-cryptographic hashes: CRC-32, Adler-32, FNV-1a.

Fast integrity checks and hash-table mixers implemented from their definitions. CRC-32
(IEEE 802.3 polynomial) is the workhorse of Ethernet, gzip, and PNG; Adler-32 is the
cheaper Zlib checksum; FNV-1a is a simple, well-dispersed hash for tables and bloom
filters. These are *not* cryptographic. Inputs are bytes or a str (UTF-8 encoded). Pure
standard library.
"""

_CRC32_POLY = 0xEDB88320
_CRC32_TABLE = None


def _to_bytes(data):
    if isinstance(data, str):
        return data.encode("utf-8")
    return bytes(data)


def _crc_table():
    global _CRC32_TABLE
    if _CRC32_TABLE is None:
        table = []
        for n in range(256):
            c = n
            for _ in range(8):
                c = (c >> 1) ^ _CRC32_POLY if (c & 1) else (c >> 1)
            table.append(c)
        _CRC32_TABLE = table
    return _CRC32_TABLE


def crc32(data):
    """CRC-32 checksum (IEEE 802.3 reflected polynomial), matching :func:`zlib.crc32`.

    Table-driven byte-at-a-time computation over the bytes of ``data`` (a ``str`` is
    UTF-8 encoded). Returns a 32-bit unsigned integer.
    """
    table = _crc_table()
    crc = 0xFFFFFFFF
    for b in _to_bytes(data):
        crc = (crc >> 8) ^ table[(crc ^ b) & 0xFF]
    return crc ^ 0xFFFFFFFF


def adler32(data):
    """Adler-32 checksum (as in zlib), matching :func:`zlib.adler32`.

    Two running sums modulo 65521 combined into a 32-bit value. Cheaper than CRC-32 but
    weaker on short inputs. Returns a 32-bit unsigned integer.
    """
    MOD = 65521
    a = 1
    b = 0
    for byte in _to_bytes(data):
        a = (a + byte) % MOD
        b = (b + a) % MOD
    return (b << 16) | a


def fnv1a_32(data):
    """32-bit FNV-1a hash: XOR then multiply per byte.

    A fast, well-dispersed non-cryptographic hash for tables and bloom filters. Returns a
    32-bit unsigned integer; ``fnv1a_32("")`` is the FNV offset basis ``2166136261``.
    """
    h = 0x811C9DC5                      # FNV offset basis
    prime = 0x01000193                  # FNV prime
    for byte in _to_bytes(data):
        h ^= byte
        h = (h * prime) & 0xFFFFFFFF
    return h
