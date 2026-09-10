"""Independent pure-Python implementations of what the C++ core computes.

These exist purely to check the bindings and the algorithms against something
that shares no code with them, so a test failure points at MOODS rather than at
a shared misunderstanding. They are written for clarity, not speed.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

ALPHABET = "ACGT"


def codes(seq: str, alphabet: str = ALPHABET) -> list[int | None]:
    """Symbol indices for a sequence; ``None`` for anything outside the alphabet."""
    lookup: dict[str, int] = {}
    for i, symbol in enumerate(alphabet):
        lookup[symbol] = i
        lookup[symbol.lower()] = i
    return [lookup.get(c) for c in seq]


def flat_bg(alphabet_size: int = 4) -> list[float]:
    return [1.0 / alphabet_size] * alphabet_size


def background(seq: str, ps: float = 1.0, alphabet: str = ALPHABET) -> list[float]:
    """``bg[j] = (count[j] + ps) / (valid_total + a * ps)``."""
    a = len(alphabet)
    counts = [0] * a
    for code in codes(seq, alphabet):
        if code is not None:
            counts[code] += 1
    total = sum(counts)
    return [(counts[j] + ps) / (total + a * ps) for j in range(a)]


def log_odds(
    counts: Sequence[Sequence[float]],
    bg: Sequence[float],
    ps: float,
    log_base: float | None = None,
) -> list[list[float]]:
    """``log(p[j][i] / bg[j])`` with ``ps`` spread over each column by ``bg``."""
    a, n = len(counts), len(counts[0])
    out = [[0.0] * n for _ in range(a)]
    for i in range(n):
        column_sum = sum(counts[j][i] + ps * bg[j] for j in range(a))
        for j in range(a):
            value = math.log((counts[j][i] + ps * bg[j]) / column_sum) - math.log(bg[j])
            out[j][i] = value if log_base is None else value / math.log(log_base)
    return out


def reverse_complement(
    matrix: Sequence[Sequence[float]], alphabet_size: int = 4
) -> list[list[float]]:
    """Reverse the columns and complement the rows.

    For a high-order matrix the row index encodes a q-tuple, so complementing a
    row means complementing and reversing the symbols it packs.
    """
    rows, cols = len(matrix), len(matrix[0])
    q = q_gram_size(rows, alphabet_size)
    shift = bits_per_symbol(alphabet_size)
    out = [[0.0] * cols for _ in range(rows)]
    for j in range(rows):
        rc = 0
        for k in range(q):
            symbol = (j >> ((q - k - 1) * shift)) & ((1 << shift) - 1)
            rc |= (alphabet_size - symbol - 1) << (k * shift)
        for i in range(cols):
            out[rc][cols - i - 1] = matrix[j][i]
    return out


def bits_per_symbol(alphabet_size: int) -> int:
    """Bits needed per symbol: ``ceil(log2(alphabet_size))``."""
    shift, size = 0, 1
    while size < alphabet_size:
        shift += 1
        size <<= 1
    return shift


def q_gram_size(rows: int, alphabet_size: int) -> int:
    """Number of symbols each matrix row indexes."""
    q, size = 0, 1
    while size < rows:
        q += 1
        size *= alphabet_size
    return q


def motif_length(matrix: Sequence[Sequence[float]], alphabet_size: int = 4) -> int:
    return len(matrix[0]) + q_gram_size(len(matrix), alphabet_size) - 1


def score_at(
    seq: str,
    pos: int,
    matrix: Sequence[Sequence[float]],
    alphabet: str = ALPHABET,
) -> float | None:
    """Score the window starting at ``pos``; ``None`` if it runs off the end or
    contains a symbol outside the alphabet."""
    a = len(alphabet)
    rows, cols = len(matrix), len(matrix[0])
    q = q_gram_size(rows, a)
    shift = bits_per_symbol(a)
    length = cols + q - 1
    if pos < 0 or pos + length > len(seq):
        return None

    window = codes(seq[pos : pos + length], alphabet)
    if any(c is None for c in window):
        return None

    score = 0.0
    for j in range(cols):
        code = 0
        for k in range(q):
            code = (code << shift) | window[j + k]
        score += matrix[code][j]
    return score


def brute_force_scan(
    seq: str,
    matrix: Sequence[Sequence[float]],
    threshold: float,
    alphabet: str = ALPHABET,
) -> list[tuple[int, float]]:
    """Every start position whose window scores at least ``threshold``."""
    length = motif_length(matrix, len(alphabet))
    hits = []
    for pos in range(len(seq) - length + 1):
        score = score_at(seq, pos, matrix, alphabet)
        if score is not None and score >= threshold:
            hits.append((pos, score))
    return hits


def max_score(matrix: Sequence[Sequence[float]], alphabet_size: int = 4) -> float:
    """Best achievable score, by exhaustive search over symbol sequences.

    Only viable for short matrices; the C++ version uses a dynamic program.
    """
    return _extreme_score(matrix, alphabet_size, max)


def min_score(matrix: Sequence[Sequence[float]], alphabet_size: int = 4) -> float:
    return _extreme_score(matrix, alphabet_size, min)


def _extreme_score(matrix, alphabet_size, pick):
    rows, cols = len(matrix), len(matrix[0])
    q = q_gram_size(rows, alphabet_size)
    shift = bits_per_symbol(alphabet_size)
    length = cols + q - 1

    best = None
    for n in range(alphabet_size**length):
        symbols = []
        rest = n
        for _ in range(length):
            symbols.append(rest % alphabet_size)
            rest //= alphabet_size
        score = 0.0
        for j in range(cols):
            code = 0
            for k in range(q):
                code = (code << shift) | symbols[j + k]
            score += matrix[code][j]
        best = score if best is None else pick(best, score)
    return best


def assert_matrix_close(actual, expected, rel: float = 1e-9, abs_: float = 1e-12) -> None:
    """Compare two matrices row by row.

    ``pytest.approx`` does not descend into nested sequences, so comparing a
    matrix directly against it raises rather than failing usefully.
    """
    import pytest

    assert len(actual) == len(expected), (
        f"row count differs: {len(actual)} != {len(expected)}"
    )
    for i, (row_a, row_e) in enumerate(zip(actual, expected)):
        assert list(row_a) == pytest.approx(list(row_e), rel=rel, abs=abs_), (
            f"row {i} differs"
        )


def assert_agrees_ignoring_ties(found, expected, threshold, tol=1e-9):
    """Compare scanner output with a brute-force scan, ignoring borderline hits.

    The scanner adds a motif's columns up window-first and then in lookahead
    order, while :func:`brute_force_scan` goes left to right. The two totals can
    differ in the last bit, so a score sitting exactly on the threshold may land
    on either side of it. Positions within ``tol`` of the threshold are
    therefore excluded from the comparison; everything else must match exactly.
    """
    found_scores = {m.pos: m.score for m in found}
    expected_scores = dict(expected)

    def decided(scores):
        return {p for p, s in scores.items() if abs(s - threshold) > tol}

    positions = decided(found_scores) | decided(expected_scores)
    for pos in sorted(positions):
        in_found = pos in found_scores
        in_expected = pos in expected_scores
        assert in_found == in_expected, (
            f"position {pos}: scanner={'hit' if in_found else 'miss'}, "
            f"brute force={'hit' if in_expected else 'miss'}"
        )
        assert abs(found_scores[pos] - expected_scores[pos]) < 1e-9, (
            f"position {pos}: {found_scores[pos]} vs {expected_scores[pos]}"
        )
