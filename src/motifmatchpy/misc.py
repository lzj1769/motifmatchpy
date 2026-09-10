"""Bit-level helpers shared by the scanning code (mirrors ``MOODS.misc``).

Matrices for models of order > 0 have one row per *q*-gram rather than per
symbol, and the scanners address those rows by packing ``q`` symbol codes into
one integer, ``shift(a)`` bits each. These are the primitives for that packing.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from . import _alphabet
from ._util import as_sequence

__all__ = ["shift", "mask", "q_gram_size", "rc_tuple", "preprocess_seq"]


def shift(a: int) -> int:
    """Bits needed to hold one symbol of an alphabet of size ``a``.

    That is ``ceil(log2(a))``: 2 for DNA.
    """
    s, b = 0, 1
    while b < a:
        s += 1
        b <<= 1
    return s


def mask(a: int) -> int:
    """Bit mask covering one symbol of an alphabet of size ``a``."""
    b = 1
    while b < a:
        b <<= 1
    return b - 1


def q_gram_size(rows: int, a: int) -> int:
    """How many symbols each row of an ``rows``-row matrix stands for.

    A 0-order matrix has ``a`` rows and so ``q == 1``; a first-order one has
    ``a**2`` rows and ``q == 2``.
    """
    q, size = 0, 1
    while size < rows:
        q += 1
        size *= a
    return q


def rc_tuple(code: int, a: int, q: int) -> int:
    """Reverse complement of a packed ``q``-gram.

    Complementing a symbol is ``a - 1 - symbol``, and the symbols come back in
    the opposite order.
    """
    s = shift(a)
    symbol_mask = (1 << s) - 1
    rc = 0
    for i in range(q):
        symbol = (code >> ((q - i - 1) * s)) & symbol_mask
        rc |= (a - symbol - 1) << (i * s)
    return rc


def preprocess_seq(
    s: str | bytes, a: int, alphabet_map: Sequence[int] | np.ndarray
) -> list[int]:
    """Bounds of the maximal runs of ``s`` made only of alphabet symbols.

    Returns a flat list of alternating start and end offsets, so
    ``[2, 7, 10, 12]`` describes the regions ``[2, 7)`` and ``[10, 12)``.
    Anything mapping to ``a`` or above -- ``N``, say -- breaks a region, because
    a match may not span a base whose identity is unknown.
    """
    table = np.asarray(alphabet_map, dtype=np.uint8)
    codes = _alphabet.encode(as_sequence(s), table)
    return [int(x) for x in region_bounds(codes, a)]


def region_bounds(codes: np.ndarray, a: int) -> np.ndarray:
    """:func:`preprocess_seq` on already-encoded symbol codes."""
    if codes.size == 0:
        return np.empty(0, dtype=np.int64)
    valid = codes < a
    # A region starts where valid turns on and ends where it turns off, so the
    # boundaries are exactly the transitions of the padded indicator.
    padded = np.concatenate(([False], valid, [False]))
    edges = np.flatnonzero(padded[1:] != padded[:-1])
    return edges.astype(np.int64)
