"""The high-level interface: motifs with names, and hits with strands.

The modules :mod:`motifmatchpy.tools`, :mod:`motifmatchpy.scan` and
:mod:`motifmatchpy.parsers` mirror the MOODS C++ API and work in bare matrices.
This module wraps them in the objects an analysis usually wants: a
:class:`Motif` that remembers its name and order, and a :class:`MotifScanner`
that handles both strands, thresholds and reverse complements for you.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from . import io as _io
from . import parsers, tools
from . import scan as _scan
from ._types import Variant

__all__ = ["Motif", "Hit", "MotifScanner", "read_motifs", "read_motif_file"]

Matrix = Sequence[Sequence[float]]

#: Files with this suffix are read as adjacent dinucleotide models.
ADM_SUFFIX = ".adm"


def _order_of(rows: int, alphabet_size: int) -> int:
    """Model order implied by a matrix with ``rows`` rows.

    A 0-order matrix has one row per symbol, a first-order one has ``a**2``.
    """
    order, n = 0, alphabet_size
    while n < rows:
        n *= alphabet_size
        order += 1
    if n != rows:
        raise ValueError(
            f"a matrix over an alphabet of size {alphabet_size} must have "
            f"{alphabet_size}**k rows, got {rows}"
        )
    return order


@dataclass(frozen=True, slots=True)
class Motif:
    """A named scoring matrix.

    ``matrix`` holds scores (usually log-odds), one row per symbol for a
    0-order motif and one row per symbol *k*-tuple for a high-order one.
    """

    name: str
    matrix: tuple[tuple[float, ...], ...]
    alphabet_size: int = 4
    strand: str = "+"

    def __post_init__(self) -> None:
        rows = tuple(tuple(float(x) for x in row) for row in self.matrix)
        if not rows:
            raise ValueError(f"motif {self.name!r} has no rows")
        width = len(rows[0])
        if width == 0:
            raise ValueError(f"motif {self.name!r} has no columns")
        if any(len(row) != width for row in rows):
            raise ValueError(f"motif {self.name!r} has rows of differing lengths")
        if self.strand not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {self.strand!r}")
        _order_of(len(rows), self.alphabet_size)  # validates the row count
        object.__setattr__(self, "matrix", rows)

    # -- shape ---------------------------------------------------------
    @property
    def order(self) -> int:
        """0 for a plain PWM, 1 for an adjacent dinucleotide model."""
        return _order_of(len(self.matrix), self.alphabet_size)

    @property
    def width(self) -> int:
        """Number of matrix columns."""
        return len(self.matrix[0])

    @property
    def length(self) -> int:
        """Number of sequence positions the motif spans."""
        return self.width + self.order

    def __len__(self) -> int:
        return self.length

    # -- scores --------------------------------------------------------
    @property
    def max_score(self) -> float:
        """Highest score this motif can assign."""
        return tools.max_score(self.matrix, self.alphabet_size)

    @property
    def min_score(self) -> float:
        """Lowest score this motif can assign."""
        return tools.min_score(self.matrix, self.alphabet_size)

    def threshold_from_p(
        self,
        p: float,
        bg: Sequence[float] | None = None,
        precision: float | None = None,
    ) -> float:
        """Score threshold with a false-positive rate of at most ``p``."""
        if bg is None:
            bg = tools.flat_bg(self.alphabet_size)
        return tools.threshold_from_p(
            self.matrix, bg, p, a=self.alphabet_size, precision=precision
        )

    # -- transformations -----------------------------------------------
    def reverse_complement(self) -> Motif:
        """The same motif scoring the opposite strand."""
        return Motif(
            name=self.name,
            matrix=tools.reverse_complement(self.matrix, self.alphabet_size),
            alphabet_size=self.alphabet_size,
            strand="-" if self.strand == "+" else "+",
        )

    def __repr__(self) -> str:
        return (
            f"Motif(name={self.name!r}, length={self.length}, "
            f"order={self.order}, strand={self.strand!r})"
        )

    # -- construction --------------------------------------------------
    @classmethod
    def from_file(
        cls,
        filename: str | os.PathLike[str],
        *,
        name: str | None = None,
        bg: Sequence[float] | None = None,
        pseudocount: float = 0.01,
        log_base: float | None = None,
        log_odds: bool = True,
        alphabet_size: int = 4,
    ) -> Motif:
        """Read a single ``.pfm`` or ``.adm`` file.

        The format is taken from the file's suffix where possible and otherwise
        detected by trying both. With ``log_odds=False`` the file is read as a
        scoring matrix and used directly instead of being converted.
        """
        path = Path(filename)
        if bg is None:
            bg = tools.flat_bg(alphabet_size)
        matrix = _read_matrix(
            path,
            bg=bg,
            pseudocount=pseudocount,
            log_base=log_base,
            log_odds=log_odds,
            alphabet_size=alphabet_size,
        )
        return cls(
            name=name if name is not None else path.name,
            matrix=matrix,
            alphabet_size=alphabet_size,
        )


def _read_matrix(
    path: Path,
    *,
    bg: Sequence[float],
    pseudocount: float,
    log_base: float | None,
    log_odds: bool,
    alphabet_size: int,
) -> list[list[float]]:
    """Read a matrix file, trying the PFM and ADM formats as appropriate."""
    errors: list[str] = []

    def try_pfm() -> list[list[float]] | None:
        try:
            mat = (
                parsers.pfm_to_log_odds(path, bg, pseudocount, log_base)
                if log_odds
                else parsers.pfm(path)
            )
        except (RuntimeError, ValueError) as exc:
            errors.append(f"as pfm: {exc}")
            return None
        # A .adm read as a .pfm parses fine but has the wrong number of rows.
        if len(mat) != alphabet_size:
            errors.append(f"as pfm: expected {alphabet_size} rows, got {len(mat)}")
            return None
        return mat

    def try_adm() -> list[list[float]] | None:
        try:
            mat = (
                parsers.adm_to_log_odds(path, bg, pseudocount, alphabet_size, log_base)
                if log_odds
                else parsers.adm_1o_terms(path, alphabet_size)
            )
        except (RuntimeError, ValueError) as exc:
            errors.append(f"as adm: {exc}")
            return None
        if len(mat) != alphabet_size**2:
            errors.append(f"as adm: expected {alphabet_size**2} rows, got {len(mat)}")
            return None
        return mat

    order = (try_adm, try_pfm) if path.suffix.lower() == ADM_SUFFIX else (try_pfm, try_adm)
    for attempt in order:
        mat = attempt()
        if mat is not None:
            return mat
    raise ValueError(f"could not parse matrix file {path} ({'; '.join(errors)})")


def read_motif_file(filename: str | os.PathLike[str], **kwargs) -> Motif:
    """Alias for :meth:`Motif.from_file`."""
    return Motif.from_file(filename, **kwargs)


def read_motifs(
    filenames: Iterable[str | os.PathLike[str]],
    *,
    bg: Sequence[float] | None = None,
    pseudocount: float = 0.01,
    log_base: float | None = None,
    log_odds: bool = True,
    alphabet_size: int = 4,
    names: Sequence[str] | None = None,
) -> list[Motif]:
    """Read several matrix files into :class:`Motif` objects."""
    filenames = list(filenames)
    if names is not None and len(names) != len(filenames):
        raise ValueError(
            f"got {len(filenames)} files but {len(names)} names"
        )
    return [
        Motif.from_file(
            filename,
            name=None if names is None else names[i],
            bg=bg,
            pseudocount=pseudocount,
            log_base=log_base,
            log_odds=log_odds,
            alphabet_size=alphabet_size,
        )
        for i, filename in enumerate(filenames)
    ]


@dataclass(frozen=True, slots=True)
class Hit:
    """One motif occurrence.

    ``variants`` is empty for an ordinary hit. It is populated for a hit that
    exists only once one or more sequence variants are applied, in which case
    :attr:`pos` is relative to the *modified* sequence.
    """

    motif: Motif
    pos: int
    strand: str
    score: float
    sequence_name: str = ""
    variants: tuple[Variant, ...] = field(default=())

    @property
    def name(self) -> str:
        """Name of the motif that matched."""
        return self.motif.name

    @property
    def end(self) -> int:
        """Position one past the last matched base."""
        return self.pos + self.motif.length

    def matched_sequence(self, seq: str) -> str:
        """The stretch of ``seq`` this hit covers, as it appears in ``seq``."""
        return seq[self.pos : self.end]

    def variant_sequence(self, seq: str) -> str:
        """The matched stretch with this hit's variant alleles applied.

        Substituted bases are upper case against a lower-case background. Only
        substitutions are rendered: an insertion or deletion shifts every
        following position, so the modified sequence no longer lines up with
        ``seq`` and cannot be reconstructed from a single window of it.
        """
        chars = list(self.matched_sequence(seq).lower())
        for variant in self.variants:
            offset = variant.start_pos - self.pos
            is_substitution = (
                variant.end_pos == variant.start_pos + 1
                and len(variant.modified_seq) == 1
            )
            if is_substitution and 0 <= offset < len(chars):
                chars[offset] = variant.modified_seq.upper()
        return "".join(chars)

    def __repr__(self) -> str:
        seq = f", sequence_name={self.sequence_name!r}" if self.sequence_name else ""
        var = f", variants={len(self.variants)}" if self.variants else ""
        return (
            f"Hit(name={self.motif.name!r}, pos={self.pos}, "
            f"strand={self.strand!r}, score={self.score:.6g}{seq}{var})"
        )


class MotifScanner:
    """Scans sequences for a fixed set of motifs.

    Preprocessing the matrices dominates the cost of a small scan, so build one
    scanner and reuse it across sequences.

    Exactly one of ``p_value`` and ``threshold`` sets the score cutoff; pass
    neither if you only intend to call :meth:`scan_best_hits`.

    ``bg`` is the background distribution used to turn a p-value into a score
    threshold. It may be a sequence of probabilities, ``"flat"`` for a uniform
    distribution, or ``"auto"`` to re-estimate it from each sequence -- which
    means recomputing thresholds per sequence, so it is markedly slower.
    """

    def __init__(
        self,
        motifs: Sequence[Motif],
        *,
        p_value: float | None = None,
        threshold: float | None = None,
        bg: Sequence[float] | str = "flat",
        both_strands: bool = True,
        window_size: int = 7,
        threshold_precision: float | None = None,
    ) -> None:
        motifs = list(motifs)
        if not motifs:
            raise ValueError("need at least one motif")
        if p_value is not None and threshold is not None:
            raise ValueError("give either p_value or threshold, not both")

        alphabet_sizes = {m.alphabet_size for m in motifs}
        if len(alphabet_sizes) != 1:
            raise ValueError(
                f"all motifs must share an alphabet size, got {sorted(alphabet_sizes)}"
            )
        self._alphabet_size = alphabet_sizes.pop()

        if isinstance(bg, str):
            if bg == "flat":
                self._bg: list[float] | None = tools.flat_bg(self._alphabet_size)
            elif bg == "auto":
                if p_value is None:
                    raise ValueError("bg='auto' only makes sense together with p_value")
                self._bg = None
            else:
                raise ValueError(f"bg must be a sequence, 'flat' or 'auto', got {bg!r}")
        else:
            self._bg = [float(x) for x in bg]

        self._motifs = tuple(motifs)
        self._both_strands = bool(both_strands)
        self._window_size = int(window_size)
        self._p_value = p_value
        self._threshold = threshold
        self._precision = threshold_precision

        # Forward motifs first, then reverse complements, so results line up
        # with self._motifs by index modulo len(self._motifs).
        self._all_motifs: tuple[Motif, ...] = self._motifs
        if self._both_strands:
            self._all_motifs += tuple(m.reverse_complement() for m in self._motifs)
        self._matrices = [m.matrix for m in self._all_motifs]

        self._cached: tuple[_scan.Scanner, list[float]] | None = None

    # -- introspection --------------------------------------------------
    @property
    def motifs(self) -> tuple[Motif, ...]:
        """The motifs this scanner was built with (forward strand)."""
        return self._motifs

    def __len__(self) -> int:
        return len(self._motifs)

    def __repr__(self) -> str:
        cutoff = (
            f"p_value={self._p_value}"
            if self._p_value is not None
            else f"threshold={self._threshold}"
            if self._threshold is not None
            else "no cutoff"
        )
        return (
            f"MotifScanner({len(self._motifs)} motifs, {cutoff}, "
            f"both_strands={self._both_strands})"
        )

    def thresholds(self, seq: str | bytes | None = None) -> list[float]:
        """Score thresholds in use, one per motif (forward strand)."""
        _, thresholds = self._prepare(seq)
        return thresholds[: len(self._motifs)]

    # -- scanning --------------------------------------------------------
    def scan(
        self,
        seq: str | bytes,
        *,
        sequence_name: str = "",
        max_hits: int | None = None,
        include_variants: bool = False,
    ) -> list[Hit]:
        """Find every occurrence in ``seq``, sorted by position.

        ``max_hits`` abandons any motif that exceeds that many occurrences,
        which keeps a too-permissive threshold from swamping the run.

        ``include_variants`` additionally reports occurrences that exist only
        because of an IUPAC ambiguity code in the sequence; those hits carry
        the indices of the variants involved in :attr:`Hit.variants`.
        """
        seq = _as_str(seq)
        scanner, _ = self._prepare(seq)

        if max_hits is None:
            results = scanner.scan(seq)
        else:
            results = scanner.scan_max_hits(seq, max_hits)

        hits = [
            Hit(
                motif=self._motifs[i % len(self._motifs)],
                pos=match.pos,
                strand=motif.strand,
                score=match.score,
                sequence_name=sequence_name,
            )
            for i, (motif, matches) in enumerate(
                zip(self._all_motifs, results, strict=True)
            )
            for match in matches
        ]

        if include_variants:
            variants = tools.snp_variants(seq)
            if variants:
                var_results = scanner.variant_matches(seq, variants)
                hits += [
                    Hit(
                        motif=self._motifs[i % len(self._motifs)],
                        pos=match.pos,
                        strand=motif.strand,
                        score=match.score,
                        sequence_name=sequence_name,
                        variants=tuple(variants[j] for j in match.variants),
                    )
                    for i, (motif, matches) in enumerate(
                        zip(self._all_motifs, var_results, strict=True)
                    )
                    for match in matches
                ]

        hits.sort(key=lambda h: (h.pos, h.motif.name, h.strand))
        return hits

    def scan_fasta(
        self,
        filename: str | os.PathLike[str],
        **kwargs,
    ) -> Iterator[Hit]:
        """Scan every record of a FASTA (or plain-text) sequence file.

        Records are read and scanned one at a time, so this streams rather than
        holding every hit in memory at once.
        """
        kwargs.pop("sequence_name", None)
        for name, seq in _io.read_sequences(filename):
            yield from self.scan(seq, sequence_name=name, **kwargs)

    def scan_best_hits(
        self,
        seq: str | bytes,
        target: int,
        *,
        sequence_name: str = "",
        iterations: int = 10,
        mult: int = 2,
        limit_mult: int = 10,
    ) -> list[Hit]:
        """Return roughly the ``target`` highest-scoring occurrences per motif.

        No threshold is needed -- one is searched for -- so this works on a
        scanner built without ``p_value`` or ``threshold``.
        """
        seq = _as_str(seq)
        results = _scan.scan_best_hits_dna(
            seq,
            self._matrices,
            target,
            iterations,
            mult,
            limit_mult,
            self._window_size,
        )
        hits = [
            Hit(
                motif=self._motifs[i % len(self._motifs)],
                pos=match.pos,
                strand=motif.strand,
                score=match.score,
                sequence_name=sequence_name,
            )
            for i, (motif, matches) in enumerate(
                zip(self._all_motifs, results, strict=True)
            )
            for match in matches
        ]
        hits.sort(key=lambda h: (h.pos, h.motif.name, h.strand))
        return hits

    def count(self, seq: str | bytes, *, max_hits: int | None = None) -> list[int]:
        """Occurrences per motif, in the order of :attr:`motifs`.

        Counts both strands when the scanner does. Cheaper than :meth:`scan`
        because no hit objects are built.
        """
        seq = _as_str(seq)
        scanner, _ = self._prepare(seq)
        if max_hits is None:
            counts = [len(r) for r in scanner.scan(seq)]
        else:
            counts = list(scanner.counts_max_hits(seq, max_hits))
        n = len(self._motifs)
        if self._both_strands:
            return [counts[i] + counts[i + n] for i in range(n)]
        return counts[:n]

    # -- internals -------------------------------------------------------
    def _prepare(self, seq: str | bytes | None) -> tuple[_scan.Scanner, list[float]]:
        """Return a scanner and its thresholds, rebuilding only when needed."""
        if self._bg is not None and self._cached is not None:
            return self._cached

        if self._bg is not None:
            bg = self._bg
        else:
            if seq is None:
                raise ValueError(
                    "bg='auto' needs a sequence; call scan() rather than thresholds()"
                )
            bg = tools.bg_from_sequence_dna(_as_str(seq), 1.0)

        thresholds = self._compute_thresholds(bg)
        scanner = _scan.Scanner(self._window_size)
        scanner.set_motifs(self._matrices, bg, thresholds)

        if self._bg is not None:
            self._cached = (scanner, thresholds)
        return scanner, thresholds

    def _compute_thresholds(self, bg: Sequence[float]) -> list[float]:
        if self._threshold is not None:
            return [float(self._threshold)] * len(self._all_motifs)
        if self._p_value is None:
            raise ValueError(
                "this scanner was built without p_value or threshold; "
                "only scan_best_hits() is available"
            )
        return [
            tools.threshold_from_p(
                m.matrix,
                bg,
                self._p_value,
                a=self._alphabet_size,
                precision=self._precision,
            )
            for m in self._all_motifs
        ]


def _as_str(seq: str | bytes) -> str:
    from ._util import as_sequence

    return as_sequence(seq)
