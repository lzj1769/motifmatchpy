"""motifmatchpy -- position weight matrix matching for DNA.

A pure-Python implementation of the MOODS algorithms (Korhonen et al.), with
the inner loops compiled by numba. It scans DNA for occurrences of position
weight matrices, including first- and higher-order models, and can report the
matches that only appear once sequence variants are applied.

Two layers are available:

* the high-level API -- :class:`Motif`, :class:`MotifScanner`, :class:`Hit` --
  which keeps track of names, strands and thresholds for you;
* :mod:`motifmatchpy.tools`, :mod:`motifmatchpy.scan`,
  :mod:`motifmatchpy.parsers` and :mod:`motifmatchpy.misc`, which mirror the
  MOODS modules of the same names and work in bare matrices.

Getting started::

    import motifmatchpy as mm

    motifs = mm.read_motifs(["MA0001.pfm", "MA0002.pfm"])
    scanner = mm.MotifScanner(motifs, p_value=1e-4)
    for hit in scanner.scan_fasta("chr1.fa"):
        print(hit.name, hit.pos, hit.strand, hit.score)

Results agree with the original MOODS C++ implementation bit for bit; the test
suite checks that against a build of it kept in ``reference/``.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _installed_version

from . import io, misc, parsers, scan, tools
from ._types import Match, MatchWithVariant, Variant
from .api import Hit, Motif, MotifScanner, read_motif_file, read_motifs
from .io import read_fasta, read_sequences
from .scan import Scanner, naive_scan_dna, scan_best_hits_dna, scan_dna
from .tools import flat_bg, threshold_from_p

try:
    __version__ = _installed_version("motifmatchpy")
except PackageNotFoundError:  # pragma: no cover - a source tree with no install
    __version__ = "0.0.0+unknown"

#: The MOODS release these algorithms were ported from.
MOODS_VERSION = "1.9.4.1"

__all__ = [
    "__version__",
    "MOODS_VERSION",
    # high-level API
    "Motif",
    "MotifScanner",
    "Hit",
    "read_motifs",
    "read_motif_file",
    # sequence IO
    "read_fasta",
    "read_sequences",
    # MOODS-level types and functions
    "Match",
    "MatchWithVariant",
    "Variant",
    "Scanner",
    "scan_dna",
    "scan_best_hits_dna",
    "naive_scan_dna",
    "flat_bg",
    "threshold_from_p",
    # modules
    "io",
    "misc",
    "parsers",
    "scan",
    "tools",
]
