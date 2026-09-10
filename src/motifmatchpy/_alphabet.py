"""Mapping sequence characters to matrix row indices.

Every scan starts by turning a sequence into an array of symbol codes: ``0`` to
``a - 1`` for the symbols the matrices are defined over, and ``a`` itself as the
"not part of the alphabet" marker that splits the sequence into scannable
regions.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

__all__ = ["DNA", "INVALID", "build_map", "encode", "alphabet_size"]

#: The default alphabet: entry *i* lists the characters read as row *i*.
DNA: tuple[str, ...] = ("aA", "cC", "gG", "tT")

#: Placeholder for "no code assigned"; :func:`build_map` replaces it with ``a``.
INVALID = 255

_TABLE_SIZE = 256


def alphabet_size(alphabet: Sequence[str] = DNA) -> int:
    return len(alphabet)


def build_map(alphabet: Sequence[str] = DNA) -> np.ndarray:
    """A 256-entry lookup table from byte value to symbol code.

    Characters not listed in ``alphabet`` map to ``len(alphabet)``, the code the
    scanners treat as breaking a scannable region.
    """
    a = len(alphabet)
    if a == 0:
        raise ValueError("alphabet must have at least one symbol")
    if a > INVALID:
        raise ValueError(f"alphabet of size {a} is too large")

    table = np.full(_TABLE_SIZE, a, dtype=np.uint8)
    for code, symbols in enumerate(alphabet):
        for char in symbols:
            value = ord(char)
            if value >= _TABLE_SIZE:
                raise ValueError(
                    f"alphabet symbol {char!r} is not a single byte; "
                    "sequences and alphabets must be ASCII"
                )
            table[value] = code
    return table


def encode(seq: str, table: np.ndarray) -> np.ndarray:
    """Turn a sequence into an array of symbol codes.

    ``seq`` must already be ASCII -- :func:`motifmatchpy._util.as_sequence`
    enforces that at the public boundary.
    """
    if not seq:
        return np.empty(0, dtype=np.uint8)
    raw = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    return table[raw]
