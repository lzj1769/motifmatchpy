"""Matrix file parsers.

Two formats are understood, both plain whitespace-separated numbers:

``.pfm``
    A JASPAR-style count matrix: one row per symbol (A, C, G, T), one column
    per motif position.

``.adm``
    An adjacent dinucleotide model: ``a * a`` rows of first-order conditional
    terms, followed by ``a`` rows giving the zero-order terms for the first
    position.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from pathlib import Path

from . import tools
from ._util import as_path, check_bg

__all__ = [
    "pfm",
    "pfm_to_log_odds",
    "adm_1o_terms",
    "adm_0o_terms",
    "adm_to_log_odds",
    "read_table",
    "ParseError",
]


class ParseError(ValueError):
    """A matrix file could not be read in the requested format."""


def read_table(filename: str | os.PathLike[str]) -> list[list[float]]:
    """Read whitespace-separated numbers into a list of rows.

    A row stops at the first token that is not a number, and rows holding no
    numbers at all are dropped -- so blank lines and comment lines disappear
    rather than becoming empty rows. This mirrors how the C++ reader consumes a
    line with a stream of ``double`` extractions.
    """
    path = Path(as_path(filename))
    rows: list[list[float]] = []
    with path.open("r", errors="replace") as handle:
        for line in handle:
            row: list[float] = []
            for token in line.split():
                try:
                    row.append(float(token))
                except ValueError:
                    break
            if row:
                rows.append(row)
    return rows


def pfm(filename: str | os.PathLike[str]) -> list[list[float]]:
    """Read a count/frequency matrix without transforming it."""
    return _read_rectangular(filename, "pfm")


def pfm_to_log_odds(
    filename: str | os.PathLike[str],
    bg: Sequence[float],
    pseudocount: float = 0.01,
    log_base: float | None = None,
) -> list[list[float]]:
    """Read a count matrix and convert it to a log-odds scoring matrix."""
    counts = _read_rectangular(filename, "pfm")
    bg = check_bg(bg, len(counts))
    return tools.log_odds(counts, bg, pseudocount, log_base)


def adm_1o_terms(
    filename: str | os.PathLike[str], a: int = 4
) -> list[list[float]]:
    """The ``a * a`` first-order conditional terms of an ADM file."""
    return _read_adm(filename, a)[: a * a]


def adm_0o_terms(
    filename: str | os.PathLike[str], a: int = 4
) -> list[list[float]]:
    """The ``a`` zero-order terms of an ADM file, for the first position."""
    return _read_adm(filename, a)[a * a :]


def adm_to_log_odds(
    filename: str | os.PathLike[str],
    bg: Sequence[float],
    pseudocount: float = 0.01,
    a: int = 4,
    log_base: float | None = None,
) -> list[list[float]]:
    """Read an adjacent dinucleotide model and convert it to log-odds scores."""
    adm = _read_adm(filename, a)
    bg = check_bg(bg, a)
    conditional = adm[: a * a]
    # Only the leading value of each zero-order row is used: it is the marginal
    # for the very first position, which has no preceding base to condition on.
    zero_order = [[adm[a * a + j][0] for j in range(a)]]
    return tools.log_odds_high_order(
        conditional, zero_order, bg, pseudocount, a, log_base
    )


def _read_rectangular(
    filename: str | os.PathLike[str], what: str
) -> list[list[float]]:
    """Read a table and require every row to be the same length."""
    rows = read_table(filename)
    if not rows or not rows[0]:
        raise ParseError(f"could not parse {what} matrix file: {filename}")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ParseError(
            f"could not parse {what} matrix file {filename}: "
            f"rows have differing lengths"
        )
    return rows


def _read_adm(filename: str | os.PathLike[str], a: int) -> list[list[float]]:
    """Read an ADM file and check its shape.

    An ADM has ``a * a`` equal-length rows of conditional terms followed by
    ``a`` rows of zero-order terms. Only the first value of each zero-order row
    is ever read, so those rows need only be non-empty.
    """
    if a < 1:
        raise ValueError(f"alphabet size must be >= 1, got {a}")
    rows = read_table(filename)
    expected = a * a + a
    if len(rows) != expected:
        raise ParseError(
            f"could not parse adm matrix file {filename}: "
            f"expected {expected} rows for an alphabet of size {a}, got {len(rows)}"
        )
    width = len(rows[0])
    if any(len(row) != width for row in rows[: a * a]):
        raise ParseError(
            f"could not parse adm matrix file {filename}: "
            f"conditional rows have differing lengths"
        )
    if any(not row for row in rows[a * a :]):
        raise ParseError(
            f"could not parse adm matrix file {filename}: "
            f"a zero-order row is empty"
        )
    return rows
