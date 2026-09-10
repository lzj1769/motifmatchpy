"""Small shared helpers."""

from __future__ import annotations

import os
from collections.abc import Sequence
from pathlib import Path

__all__ = ["as_sequence", "as_path", "check_bg"]


def as_sequence(seq: str | bytes | bytearray | memoryview) -> str:
    """Normalise a nucleotide sequence to an ASCII ``str``.

    Sequences frequently arrive as ``bytes`` (from ``mmap``, pysam, or a binary
    file read), so accept those rather than forcing the caller to decode.

    Non-ASCII input is rejected here. Sequences are looked up in a 256-entry
    byte table to turn characters into symbol codes, and anything outside ASCII
    has no byte to look up. Checking costs nothing -- CPython already knows
    whether a string is ASCII -- and the alternative is a confusing failure
    much deeper in.
    """
    if isinstance(seq, str):
        if not seq.isascii():
            raise ValueError(
                "sequence must be ASCII; found a character outside it "
                f"at index {_first_non_ascii(seq)}"
            )
        return seq
    if isinstance(seq, (bytes, bytearray, memoryview)):
        return bytes(seq).decode("ascii")
    raise TypeError(f"sequence must be str or bytes-like, got {type(seq).__name__}")


def _first_non_ascii(seq: str) -> int:
    for i, c in enumerate(seq):
        if not c.isascii():
            return i
    return -1


def as_path(filename: str | os.PathLike[str]) -> str:
    """Normalise a path and fail early -- and clearly -- if it is not readable.

    The C++ parsers report an unreadable file the same way they report a
    malformed one, which makes a typo in a path needlessly hard to diagnose.
    """
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(f"no such file: {path}")
    if path.is_dir():
        raise IsADirectoryError(f"expected a file, got a directory: {path}")
    return str(path)


def check_bg(bg: Sequence[float], alphabet_size: int = 4) -> list[float]:
    """Validate a background distribution."""
    bg = [float(x) for x in bg]
    if len(bg) != alphabet_size:
        raise ValueError(
            f"background must have {alphabet_size} entries, got {len(bg)}"
        )
    if any(x < 0 for x in bg):
        raise ValueError(f"background probabilities must be non-negative: {bg}")
    if not bg or sum(bg) <= 0:
        raise ValueError(f"background probabilities must sum to a positive value: {bg}")
    return bg
