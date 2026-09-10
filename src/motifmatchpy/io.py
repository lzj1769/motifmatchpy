"""Reading sequence files.

Handles FASTA and plain-text sequence files, transparently gzip-compressed.
"""

from __future__ import annotations

import gzip
import io as _io
import os
from collections.abc import Iterator
from pathlib import Path

__all__ = ["read_fasta", "read_sequences", "open_text"]

_GZIP_MAGIC = b"\x1f\x8b"


def open_text(filename: str | os.PathLike[str]) -> _io.TextIOBase:
    """Open a possibly gzip-compressed text file.

    Compression is detected from the file's magic bytes rather than its name,
    so a gzipped file called ``.fa`` still works.
    """
    path = Path(filename)
    with open(path, "rb") as probe:
        compressed = probe.read(2) == _GZIP_MAGIC
    if compressed:
        return gzip.open(path, "rt")
    return open(path)


def read_fasta(
    filename: str | os.PathLike[str], *, full_header: bool = False
) -> Iterator[tuple[str, str]]:
    """Yield ``(name, sequence)`` for each record in a FASTA file.

    The name is the sequence identifier: the header up to the first whitespace,
    following the usual convention that anything after it is free-text
    description. Pass ``full_header=True`` to keep the whole header line.
    """
    with open_text(filename) as handle:
        name: str | None = None
        chunks: list[str] = []
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = _record_name(line, full_header)
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks)


def _record_name(header_line: str, full_header: bool) -> str:
    header = header_line[1:].strip()
    if full_header:
        return header
    return header.split(maxsplit=1)[0] if header else ""


def read_sequences(
    filename: str | os.PathLike[str], *, full_header: bool = False
) -> Iterator[tuple[str, str]]:
    """Yield ``(name, sequence)`` from a FASTA *or* plain-text sequence file.

    A plain-text file is treated as a single record named after the file.
    """
    path = Path(filename)
    with open_text(path) as handle:
        for line in handle:
            if line.strip():
                is_fasta = line.startswith(">")
                break
        else:
            return  # empty file: no records

    if is_fasta:
        yield from read_fasta(path, full_header=full_header)
    else:
        with open_text(path) as handle:
            yield path.name, "".join(line.strip() for line in handle)
