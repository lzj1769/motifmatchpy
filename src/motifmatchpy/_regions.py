"""BED intervals and reference sequences for genomic motif scanning."""

from __future__ import annotations

import os
import re
from collections.abc import Iterator
from dataclasses import dataclass

from .io import open_text


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int


def read_peaks(filename: str | os.PathLike[str]) -> dict[str, list[Region]]:
    """Read BED3+ (including narrowPeak), retaining duplicate/overlapping peaks."""
    regions: dict[str, list[Region]] = {}
    with open_text(filename) as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields or fields[0].startswith("#") or fields[0] in {"track", "browser"}:
                continue
            try:
                chrom, start, end = fields[:3]
                start, end = int(start), int(end)
                if start < 0 or end <= start:
                    raise ValueError
            except ValueError as exc:
                raise ValueError(
                    f"{filename}:{number}: expected BED chrom start end with 0 <= start < end"
                ) from exc
            regions.setdefault(chrom, []).append(Region(chrom, start, end))
    for peaks in regions.values():
        peaks.sort(key=lambda region: (region.start, region.end))
    return regions


def peak_sequences(
    genome: str | os.PathLike[str], regions: dict[str, list[Region]],
) -> Iterator[tuple[Region, str]]:
    """Stream reference FASTA, keeping at most one requested chromosome in memory.

    No FASTA index or native dependency is needed. Output follows FASTA contig
    order and then peak coordinates. Intervals must fit the reference exactly;
    chromosome aliases and out-of-bounds clipping are deliberately not guessed.
    """
    if not regions:
        return
    seen: set[str] = set()
    name: str | None = None
    chunks: list[str] = []

    def emit() -> Iterator[tuple[Region, str]]:
        if name not in regions:
            return
        sequence = "".join(chunks)
        for region in regions[name]:
            if region.end > len(sequence):
                raise ValueError(
                    f"peak {region.chrom}:{region.start}-{region.end} exceeds "
                    f"reference length {len(sequence)}"
                )
            yield region, sequence[region.start:region.end]

    with open_text(genome) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith(">"):
                yield from emit()
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"{genome}:{number}: empty FASTA header")
                name = fields[0]
                if name in seen:
                    raise ValueError(f"{genome}: duplicate FASTA identifier {name!r}")
                seen.add(name)
                chunks = []
            elif line.strip():
                if name is None:
                    raise ValueError(f"{genome}:{number}: expected a FASTA header")
                if name in regions:
                    sequence_line = line.strip()
                    if not sequence_line.isascii() or re.search(r"\s", sequence_line):
                        raise ValueError(
                            f"{genome}:{number}: invalid FASTA sequence whitespace/ASCII"
                        )
                    chunks.append(sequence_line)
        yield from emit()
    missing = regions.keys() - seen
    if missing:
        raise ValueError(f"chromosomes missing from reference FASTA: {', '.join(sorted(missing))}")
