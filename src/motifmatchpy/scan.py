"""Scanning sequences for motif occurrences.

:class:`Scanner` is the workhorse: give it matrices and thresholds once, then
scan as many sequences as you like. The module-level functions are convenience
wrappers that build a scanner for a single call.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from . import _alphabet, _kernels, misc, tools
from ._kernels import ALL_HITS, COUNTS, MAX_HITS, MotifTables
from ._motif import build_motif
from ._types import Match, MatchWithVariant, Variant
from ._util import as_sequence, check_bg

__all__ = [
    "Match",
    "MatchWithVariant",
    "Variant",
    "Scanner",
    "scan_dna",
    "scan",
    "scan_best_hits_dna",
    "naive_scan_dna",
]

Matrix = Sequence[Sequence[float]]

DEFAULT_WINDOW_SIZE = 7


class Scanner:
    """Scans sequences for a fixed set of motifs.

    Preparing the matrices costs far more than scanning a short sequence, so
    build one scanner and reuse it rather than calling :func:`scan_dna` in a
    loop.

    ``window_size`` trades preparation for speed: the lookup table has one entry
    per possible window content, so it grows as ``alphabet_size ** window_size``
    while longer windows filter more aggressively. The default of 7 is a good
    balance for DNA and rarely worth changing.
    """

    def __init__(
        self,
        window_size: int = DEFAULT_WINDOW_SIZE,
        alphabet: Sequence[str] | None = None,
    ) -> None:
        window_size = int(window_size)
        if window_size < 1:
            raise ValueError(f"window size must be at least 1, got {window_size}")

        self.alphabet: tuple[str, ...] = tuple(
            _alphabet.DNA if alphabet is None else alphabet
        )
        self.a = len(self.alphabet)
        self.l = window_size
        self._shift = misc.shift(self.a)
        self._map = _alphabet.build_map(self.alphabet)

        table_size = 1 << (self._shift * self.l)
        if table_size > 1 << 26:
            raise ValueError(
                f"window size {window_size} over an alphabet of size {self.a} "
                f"needs a lookup table of {table_size} entries; use a smaller window"
            )

        self._motifs: list = []
        self._tables: MotifTables | None = None

    # -- setup -----------------------------------------------------------
    def set_motifs(
        self,
        matrices: Sequence[Matrix],
        bg: Sequence[float],
        thresholds: Sequence[float],
    ) -> None:
        """Prepare a set of matrices for scanning.

        ``bg`` steers how the scanner lays out its search -- which slice of each
        motif goes in the window, and in what order the rest is examined. It
        does not enter the scores, so it changes how fast a scan runs but never
        what it finds.
        """
        matrices = [_as_matrix(m) for m in matrices]
        thresholds = [float(t) for t in thresholds]
        if len(matrices) != len(thresholds):
            raise ValueError(
                f"got {len(matrices)} matrices but {len(thresholds)} thresholds"
            )
        bg = np.array(check_bg(bg, self.a))

        self._motifs = [
            build_motif(matrix, bg, self.l, threshold, self.a)
            for matrix, threshold in zip(matrices, thresholds)
        ]
        self._tables = _pack(self._motifs, self.l, self._shift)

    def size(self) -> int:
        """How many motifs are loaded."""
        return len(self._motifs)

    def __len__(self) -> int:
        return len(self._motifs)

    def __repr__(self) -> str:
        return (
            f"Scanner(window_size={self.l}, alphabet_size={self.a}, "
            f"motifs={len(self._motifs)})"
        )

    @property
    def max_motif_size(self) -> int:
        """The longest motif loaded, in sequence positions."""
        return max((m.size for m in self._motifs), default=0)

    # -- scanning --------------------------------------------------------
    def scan(self, seq: str | bytes) -> list[list[Match]]:
        """Every match above threshold, as one list per motif, ordered by position."""
        motif, pos, score, _ = self._run(seq, ALL_HITS, 0)
        return _group(motif, pos, score, len(self._motifs))

    def scan_max_hits(self, seq: str | bytes, max_hits: int) -> list[list[Match]]:
        """Matches per motif, giving up on a motif once it has ``max_hits``.

        The matches kept are the first ``max_hits`` encountered, which is to say
        the leftmost ones -- not the highest scoring. This is a guard against a
        threshold so loose that one motif would swamp the run, not a way to pick
        out the best sites; :func:`scan_best_hits_dna` does that.
        """
        max_hits = int(max_hits)
        if max_hits < 0:
            raise ValueError(f"max_hits must be non-negative, got {max_hits}")
        motif, pos, score, _ = self._run(seq, MAX_HITS, max_hits)
        return _group(motif, pos, score, len(self._motifs))

    def counts_max_hits(self, seq: str | bytes, max_hits: int) -> list[int]:
        """How many matches each motif has, with the same per-motif cut-off.

        Cheaper than :meth:`scan_max_hits` because no match objects are built.
        """
        max_hits = int(max_hits)
        if max_hits < 0:
            raise ValueError(f"max_hits must be non-negative, got {max_hits}")
        _, _, _, counts = self._run(seq, COUNTS, max_hits)
        return [int(c) for c in counts]

    def _run(self, seq: str | bytes, mode: int, max_hits: int):
        if self._tables is None:
            empty_i = np.empty(0, dtype=np.int64)
            return empty_i, empty_i, np.empty(0, dtype=np.float64), empty_i
        codes = _alphabet.encode(as_sequence(seq), self._map)
        bounds = misc.region_bounds(codes, self.a)
        if bounds.size == 0:
            return (
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.float64),
                np.zeros(len(self._motifs), dtype=np.int64),
            )
        return _kernels.process_matches(
            codes, bounds, self._shift, self.l, self._tables, mode, max_hits
        )

    # -- variants --------------------------------------------------------
    def variant_matches(
        self,
        seq: str | bytes,
        variants: Sequence[Variant],
        max_depth: int = 0,
    ) -> list[list[MatchWithVariant]]:
        """Matches that exist only once one or more ``variants`` are applied.

        Each result is reported with the variants responsible for it. A match
        the unmodified sequence already contains is *not* reported here -- use
        :meth:`scan` for those -- and because insertions and deletions shift
        everything after them, positions are relative to the modified sequence.

        ``max_depth`` caps how many variants may be combined in one match; 0
        means no limit beyond what the motifs' length allows.
        """
        seq = as_sequence(seq)
        n_motifs = len(self._motifs)
        results: list[list[MatchWithVariant]] = [[] for _ in range(n_motifs)]
        if not variants or n_motifs == 0:
            return results

        # Sort by position, keeping track of where each variant came from so the
        # matches can be reported in the caller's own numbering.
        original = sorted(
            range(len(variants)),
            key=lambda i: (variants[i].start_pos, variants[i].end_pos, i),
        )
        ordered = [variants[i] for i in original]
        span = self.max_motif_size

        for i, variant in enumerate(ordered):
            if _is_edge_deletion(variant, len(seq)):
                continue
            # Start from just enough preceding sequence that a motif ending on
            # this variant still fits.
            if variant.start_pos >= span - 1:
                prefix_start = variant.start_pos - span + 1
                prefix_length = span - 1
            else:
                prefix_start = 0
                prefix_length = variant.start_pos
            prefix = seq[prefix_start : prefix_start + prefix_length] + variant.modified_seq
            self._variant_matches_from(
                results, [i], prefix, prefix_start, prefix_length, 1,
                seq, ordered, max_depth,
            )

        return [
            [
                MatchWithVariant(
                    hit.pos, hit.score, tuple(original[j] for j in hit.variants)
                )
                for hit in motif_results
            ]
            for motif_results in results
        ]

    def _variant_matches_from(
        self, results, active, prefix, seq_start, variant_start, depth,
        seq, variants, max_depth,
    ) -> None:
        """Scan the sequence built so far, then extend it with further variants."""
        span = self.max_motif_size
        # The last position that must be covered for the first active variant to
        # matter, and the first that must be covered for the last one to.
        first_required = variant_start + len(variants[active[0]].modified_seq) - 1
        last_required = len(prefix) - len(variants[active[-1]].modified_seq)

        # Enough of the original sequence to let a motif starting inside the
        # last variant finish.
        remaining = 0
        if len(prefix) < span + first_required:
            remaining = span - (len(prefix) - first_required)

        tail_start = variants[active[-1]].end_pos
        current = prefix + seq[tail_start : tail_start + remaining]

        for motif_index, matches in enumerate(self.scan(current)):
            motif_length = self._motifs[motif_index].size
            for hit in matches:
                # Keep only matches that actually cover every active variant;
                # anything else would be found by an ordinary scan.
                if hit.pos + motif_length >= last_required + 1 and hit.pos <= first_required:
                    results[motif_index].append(
                        MatchWithVariant(
                            hit.pos + seq_start, hit.score, tuple(active)
                        )
                    )

        if depth == max_depth:
            return

        last = variants[active[-1]]
        for nxt in range(active[-1] + 1, len(variants)):
            candidate = variants[nxt]
            if candidate.start_pos >= remaining + last.end_pos:
                break
            if not _variants_compatible(last, candidate):
                continue
            if _is_edge_deletion(candidate, len(seq)):
                continue
            extended = (
                prefix
                + seq[last.end_pos : candidate.start_pos]
                + candidate.modified_seq
            )
            self._variant_matches_from(
                results, [*active, nxt], extended, seq_start, variant_start,
                depth + 1, seq, variants, max_depth,
            )


def _variants_compatible(current: Variant, following: Variant) -> bool:
    """Whether ``following`` can be applied after ``current`` without overlapping.

    Two variants may abut, but only if at least one of them actually replaces
    something -- otherwise two insertions at the same point would have no
    defined order.
    """
    if current.end_pos < following.start_pos:
        return True
    return current.end_pos == following.start_pos and (
        current.start_pos < current.end_pos or following.start_pos < following.end_pos
    )


def _is_edge_deletion(variant: Variant, seq_length: int) -> bool:
    """A deletion running off either end of the sequence, which means nothing here."""
    if variant.modified_seq:
        return False
    return variant.start_pos == 0 or variant.end_pos == seq_length


# ------------------------------------------------------------------ packing
def _pack(motifs: list, window_size: int, shift: int) -> MotifTables:
    """Flatten the motifs and build the window lookup table."""
    table_size = 1 << (shift * window_size)

    kind = np.array([m.kind for m in motifs], dtype=np.int64)
    size = np.array([m.size for m in motifs], dtype=np.int64)
    window_pos = np.array([m.wp for m in motifs], dtype=np.int64)
    threshold = np.array([m.T for m in motifs], dtype=np.float64)

    mat_data, mat_index = _concatenate([m.mat.ravel() for m in motifs], np.float64)
    mat_cols = np.array([m.mat.shape[1] for m in motifs], dtype=np.int64)

    la_order, la_index = _concatenate(
        [
            m.lookahead_order if m.kind == 0 else np.empty(0, dtype=np.int64)
            for m in motifs
        ],
        np.int64,
    )
    la_score, _ = _concatenate(
        [
            m.lookahead_scores if m.kind == 0 else np.empty(0, dtype=np.float64)
            for m in motifs
        ],
        np.float64,
    )

    high = [m if m.kind == 1 else None for m in motifs]
    q = np.array([h.q if h else 1 for h in high], dtype=np.int64)
    m_shift = np.array([h.SHIFT if h else shift for h in high], dtype=np.int64)
    gram_mask = np.array([h.MASK if h else 0 for h in high], dtype=np.int64)
    ctx_shift = np.array([h.Q_SHIFT if h else 0 for h in high], dtype=np.int64)
    ctx_mask = np.array([h.Q_MASK if h else 0 for h in high], dtype=np.int64)
    ctx_size = np.array([h.Q_CODE_SIZE if h else 1 for h in high], dtype=np.int64)

    prefix_data, prefix_index = _concatenate(
        [h.P.ravel() if h else np.empty(0, dtype=np.float64) for h in high], np.float64
    )
    suffix_data, suffix_index = _concatenate(
        [h.S.ravel() if h else np.empty(0, dtype=np.float64) for h in high], np.float64
    )

    # Window lookup table: for each motif, the window contents it could match.
    codes_per_motif = []
    motif_per_entry = []
    scores_per_entry = []
    for k, motif in enumerate(motifs):
        candidate, score = motif.window_scores(shift)
        hits = np.flatnonzero(candidate)
        codes_per_motif.append(hits)
        motif_per_entry.append(np.full(hits.size, k, dtype=np.int64))
        scores_per_entry.append(score[hits])

    all_codes = _concat_or_empty(codes_per_motif, np.int64)
    all_motifs = _concat_or_empty(motif_per_entry, np.int64)
    all_scores = _concat_or_empty(scores_per_entry, np.float64)

    # Group by window code, motifs in index order within each code -- the same
    # order the C++ table is built in.
    order = np.lexsort((all_motifs, all_codes))
    hit_motif = all_motifs[order]
    hit_score = all_scores[order]
    hit_index = np.zeros(table_size + 1, dtype=np.int64)
    hit_index[1:] = np.cumsum(np.bincount(all_codes, minlength=table_size))
    hit_full = (size[hit_motif] <= window_size).astype(np.uint8)

    return MotifTables(
        hit_index=hit_index,
        hit_motif=hit_motif,
        hit_score=hit_score,
        hit_full=hit_full,
        kind=kind,
        size=size,
        window_pos=window_pos,
        threshold=threshold,
        mat_index=mat_index,
        mat_data=mat_data,
        mat_cols=mat_cols,
        la_index=la_index,
        la_order=la_order,
        la_score=la_score,
        q=q,
        shift=m_shift,
        gram_mask=gram_mask,
        ctx_shift=ctx_shift,
        ctx_mask=ctx_mask,
        ctx_size=ctx_size,
        prefix_index=prefix_index,
        prefix_data=prefix_data,
        suffix_index=suffix_index,
        suffix_data=suffix_data,
    )


def _concatenate(parts: list[np.ndarray], dtype) -> tuple[np.ndarray, np.ndarray]:
    """Lay arrays end to end and return them with their start offsets."""
    index = np.zeros(len(parts) + 1, dtype=np.int64)
    index[1:] = np.cumsum([p.size for p in parts])
    return _concat_or_empty(parts, dtype), index


def _concat_or_empty(parts: list[np.ndarray], dtype) -> np.ndarray:
    if not parts:
        return np.empty(0, dtype=dtype)
    return np.concatenate(parts).astype(dtype, copy=False)


def _group(
    motif: np.ndarray, pos: np.ndarray, score: np.ndarray, n_motifs: int
) -> list[list[Match]]:
    """Split flat match arrays into one list per motif, keeping scan order.

    A stable sort by motif leaves each motif's matches in the order the scan
    found them, which is left to right. Converting the positions and scores in
    bulk with ``tolist()`` is far cheaper than unboxing them one at a time, and
    on a permissive threshold this loop dominates the whole scan.
    """
    if motif.size == 0:
        return [[] for _ in range(n_motifs)]

    order = np.argsort(motif, kind="stable")
    by_motif = motif[order]
    positions = pos[order].tolist()
    scores = score[order].tolist()
    edges = np.searchsorted(by_motif, np.arange(n_motifs + 1))
    return [
        list(map(Match, positions[start:stop], scores[start:stop]))
        for start, stop in zip(edges[:-1], edges[1:])
    ]


# ------------------------------------------------------- one-shot functions
def scan_dna(
    seq: str | bytes,
    matrices: Sequence[Matrix],
    bg: Sequence[float],
    thresholds: Sequence[float],
    window_size: int = DEFAULT_WINDOW_SIZE,
) -> list[list[Match]]:
    """Scan a DNA sequence, returning one list of matches per matrix."""
    scanner = Scanner(window_size)
    scanner.set_motifs(matrices, bg, thresholds)
    return scanner.scan(seq)


def scan(
    seq: str | bytes,
    matrices: Sequence[Matrix],
    bg: Sequence[float],
    thresholds: Sequence[float],
    window_size: int,
    alphabet: Sequence[str],
) -> list[list[Match]]:
    """Scan a sequence over a custom alphabet.

    ``alphabet`` maps matrix rows to sequence characters: entry *i* lists every
    character that should be read as row *i*, e.g. ``['aA', 'cC', 'gG', 'tT']``.
    """
    scanner = Scanner(window_size, alphabet)
    scanner.set_motifs(matrices, bg, thresholds)
    return scanner.scan(seq)


def scan_best_hits_dna(
    seq: str | bytes,
    matrices: Sequence[Matrix],
    target: int,
    iterations: int = 10,
    MULT: int = 2,  # MOODS parameter name
    LIMIT_MULT: int = 10,  # MOODS parameter name
    window_size: int = DEFAULT_WINDOW_SIZE,
) -> list[list[Match]]:
    """Return roughly the ``target`` best-scoring matches for each matrix.

    No threshold is needed: one is searched for. Each matrix is first tried at
    the highest score it can reach; those that are still too common there are
    given up on entirely, and the rest get a threshold found by binary search,
    aiming for between ``target`` and ``MULT * target`` matches and never
    returning more than ``LIMIT_MULT * target``. ``iterations`` caps the search
    (0 for no cap).
    """
    seq = as_sequence(seq)
    target = int(target)
    if target < 1:
        raise ValueError(f"target must be >= 1, got {target}")
    matrices = [_as_matrix(m) for m in matrices]
    if not matrices:
        return []

    bg = tools.bg_from_sequence_dna(seq, 0.01)
    mult = min(int(MULT), int(LIMIT_MULT))
    limit = int(LIMIT_MULT) * target
    # An initial guess aimed at the middle of the acceptable band.
    guess_p = ((1 + int(MULT)) / 3.0) * target / max(len(seq), 1)

    n = len(matrices)
    upper = [tools.max_score(m, 4) - tools.min_delta(m) / 2 for m in matrices]
    lower = [0.0] * n
    current = list(upper)
    settled = [False] * n
    remaining = n

    def count_at(indices, thresholds):
        scanner = Scanner(window_size)
        scanner.set_motifs([matrices[i] for i in indices], bg, thresholds)
        return scanner.counts_max_hits(seq, limit)

    # Start from the best score each matrix can reach: if even that is common,
    # no threshold will do, and if it already lands in range we are done.
    counts = count_at(range(n), upper)
    for i, hits in enumerate(counts):
        if target <= hits < limit:
            settled[i] = True
            remaining -= 1
        elif hits >= limit:
            current[i] = tools.max_score(matrices[i], 4) + 1.0
            settled[i] = True
            remaining -= 1

    for i in range(n):
        if not settled[i]:
            current[i] = tools.threshold_from_p(matrices[i], bg, guess_p, a=4)
            lower[i] = tools.min_score(matrices[i], 4) - 1.0

    iteration = 0
    while remaining > 0:
        iteration += 1
        pending = [i for i in range(n) if not settled[i]]
        counts = count_at(pending, [current[i] for i in pending])

        for i, hits in zip(pending, counts):
            threshold = current[i]
            if target <= hits < mult * target:
                settled[i] = True
                remaining -= 1
            elif hits >= mult * target:
                current[i] = (upper[i] + threshold) / 2.0
                lower[i] = threshold
            else:
                current[i] = (lower[i] + threshold) / 2.0
                upper[i] = threshold

            # Give up if the search has stopped moving or run out of rounds.
            if not settled[i] and (current[i] == threshold or iteration == iterations):
                settled[i] = True
                if hits >= limit:
                    current[i] = upper[i]
                remaining -= 1

    scanner = Scanner(window_size)
    scanner.set_motifs(matrices, bg, current)
    return scanner.scan(seq)


def naive_scan_dna(
    seq: str | bytes,
    matrix: Matrix,
    threshold: float,
    a: int | None = None,
) -> list[Match]:
    """Brute-force scan of a single matrix: score every position, keep the rest.

    Slow, and useful mainly as a reference the fast scanners can be checked
    against.
    """
    matrix = np.asarray(_as_matrix(matrix), dtype=np.float64)
    threshold = float(threshold)
    alphabet_size = 4 if a is None else int(a)

    codes = _alphabet.encode(as_sequence(seq), _alphabet.build_map(_alphabet.DNA))
    bounds = misc.region_bounds(codes, 4)

    rows, cols = matrix.shape
    q = 1 if a is None else misc.q_gram_size(rows, alphabet_size)
    shift = misc.shift(alphabet_size)
    length = cols + q - 1

    results: list[Match] = []
    for r in range(0, bounds.size, 2):
        start, end = int(bounds[r]), int(bounds[r + 1])
        for i in range(start, end - length + 1):
            score = 0.0
            for j in range(cols):
                code = 0
                for k in range(q):
                    code = (code << shift) | int(codes[i + j + k])
                score += matrix[code, j]
            if score >= threshold:
                results.append(Match(i, score))
    return results


def _as_matrix(mat: Matrix) -> list[list[float]]:
    if isinstance(mat, np.ndarray):
        if mat.ndim != 2 or mat.size == 0:
            raise ValueError("matrix must be a non-empty 2-D array")
        return mat.astype(np.float64, copy=False).tolist()
    rows = [[float(x) for x in row] for row in mat]
    if not rows or not rows[0]:
        raise ValueError("matrix has no rows or no columns")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("matrix rows have differing lengths")
    return rows
