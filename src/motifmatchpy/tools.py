"""Matrix transformations, backgrounds and score thresholds.

A pure-Python port of ``core/moods_tools.cpp``. The dynamic programs that turn
a p-value into a score threshold quantise scores to integers exactly the way
MOODS does, so the two implementations return the same threshold rather than
merely a close one.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

from . import misc
from ._alphabet import DNA, build_map, encode
from ._types import Variant
from ._util import as_sequence, check_bg

__all__ = [
    "DEFAULT_DP_PRECISION",
    "Variant",
    "flat_bg",
    "bg_from_sequence_dna",
    "snp_variants",
    "reverse_complement",
    "log_odds",
    "log_odds_high_order",
    "threshold_from_p",
    "threshold_from_p_with_precision",
    "max_score",
    "min_score",
    "min_delta",
]

#: Default granularity of the dynamic program behind :func:`threshold_from_p`.
#: Scores are rounded to multiples of ``1 / DEFAULT_DP_PRECISION``.
DEFAULT_DP_PRECISION: float = 2000.0

Matrix = Sequence[Sequence[float]]

# IUPAC ambiguity codes and the bases each one stands for.
IUPAC_ALTERNATIVES: dict[str, str] = {
    "W": "AT", "S": "CG", "M": "AC", "K": "GT", "R": "AG", "Y": "CT",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
}


# --------------------------------------------------------------- backgrounds
def flat_bg(alphabet_size: int = 4) -> list[float]:
    """Uniform background distribution over ``alphabet_size`` symbols."""
    if alphabet_size < 1:
        raise ValueError(f"alphabet_size must be >= 1, got {alphabet_size}")
    return [1.0 / alphabet_size] * alphabet_size


def bg_from_sequence_dna(seq: str | bytes, ps: float = 1.0) -> list[float]:
    """Estimate an ACGT background distribution from a DNA sequence.

    ``ps`` is a pseudocount added to each of the four base counts. Anything that
    is not an unambiguous base is ignored entirely -- it contributes to neither
    the counts nor the total.
    """
    codes = encode(as_sequence(seq), build_map(DNA))
    counts = np.bincount(codes, minlength=5)[:4].astype(np.float64)
    total = counts.sum()
    ps = float(ps)
    return [float((counts[j] + ps) / (total + 4 * ps)) for j in range(4)]


def snp_variants(seq: str | bytes) -> list[Variant]:
    """One variant per alternative base of every IUPAC ambiguity code in ``seq``.

    The result is ordered by position, and within a position by the alternative
    bases in alphabetical order -- the order the scanners expect.
    """
    seq = as_sequence(seq)
    return [
        Variant(i, i + 1, base)
        for i, char in enumerate(seq)
        for base in IUPAC_ALTERNATIVES.get(char.upper(), "")
    ]


# ---------------------------------------------------------------- transforms
def reverse_complement(mat: Matrix, a: int | None = None) -> list[list[float]]:
    """Reverse complement of a score matrix.

    Reverses the columns and complements the rows. Pass ``a`` -- the alphabet
    size -- for a high-order matrix, whose rows stand for q-grams rather than
    single symbols; without it the matrix is treated as 0-order.
    """
    matrix = _as_array(mat)
    rows = matrix.shape[0]
    if a is None:
        return _to_lists(matrix[::-1, ::-1])

    q = misc.q_gram_size(rows, a)
    out = np.empty_like(matrix)
    rc_rows = np.array([misc.rc_tuple(j, a, q) for j in range(rows)])
    out[rc_rows, :] = matrix[:, ::-1]
    return _to_lists(out)


def log_odds(
    mat: Matrix,
    bg: Sequence[float],
    ps: float,
    log_base: float | None = None,
    a: int | None = None,
) -> list[list[float]]:
    """Turn a count or frequency matrix into a log-odds scoring matrix.

    Each column is normalised to a distribution after adding ``ps`` spread over
    the symbols in proportion to ``bg``, then divided by ``bg`` and logged.
    ``log_base`` defaults to the natural logarithm.

    MOODS overloads this name, spelling the high-order transform
    ``log_odds(mat, low_order_terms, bg, ps, a)``. That call is recognised here
    and forwarded to :func:`log_odds_high_order`; to combine it with a log base,
    call that function directly.
    """
    if _is_matrix(bg):
        low_order_terms, high_bg, high_ps, high_a = bg, ps, log_base, a
        if high_ps is None or high_a is None:
            raise TypeError(
                "the high-order form is log_odds(mat, low_order_terms, bg, ps, a); "
                "use log_odds_high_order() to also pass a log base"
            )
        return log_odds_high_order(mat, low_order_terms, high_bg, high_ps, high_a)

    counts = _as_array(mat)
    rows = counts.shape[0]
    bg_arr = np.array(check_bg(bg, rows if a is None else a))
    ps = float(ps)

    adjusted = counts + ps * bg_arr[:, None]
    scores = _log(adjusted / adjusted.sum(axis=0)) - _log(bg_arr)[:, None]
    return _to_lists(_rescale(scores, log_base))


def log_odds_high_order(
    mat: Matrix,
    low_order_terms: Matrix,
    bg: Sequence[float],
    ps: float,
    a: int = 4,
    log_base: float | None = None,
) -> list[list[float]]:
    """Log-odds transform of a high-order matrix.

    ``mat`` holds the conditional terms: row ``(context << shift) | symbol``
    scores that symbol given the preceding context. Each context is normalised
    independently.

    The first columns of the motif have no full context, so ``low_order_terms``
    supplies the shorter conditionals for them; those are folded into column 0,
    which is why a high-order matrix spans ``columns + q - 1`` positions.
    """
    conditional = _as_array(mat)
    rows, cols = conditional.shape
    low = _as_array(low_order_terms)
    bg_arr = np.array(check_bg(bg, a))
    ps = float(ps)

    q = misc.q_gram_size(rows, a)
    shift = misc.shift(a)
    context_count = 1 << (shift * (q - 1))

    # Highest-order terms: normalise the `a` rows sharing each context.
    grouped = conditional.reshape(context_count, a, cols) + ps * bg_arr[None, :, None]
    scores = _log(grouped / grouped.sum(axis=1, keepdims=True)) - _log(bg_arr)[
        None, :, None
    ]
    out = scores.reshape(rows, cols)

    # Lower-order terms for the leading positions, each added to column 0 of
    # every row whose q-gram starts with the corresponding prefix.
    for r in range(q - 1):
        if r >= low.shape[0]:
            raise ValueError(
                f"low_order_terms needs {q - 1} rows for a {q}-gram matrix, "
                f"got {low.shape[0]}"
            )
        prefix_count = 1 << (shift * r)
        terms = low[r][: prefix_count * a].reshape(prefix_count, a) + ps * bg_arr
        lo = _log(terms / terms.sum(axis=1, keepdims=True)) - _log(bg_arr)
        # Each (prefix, symbol) pair fixes the top bits of the row index; every
        # row sharing those bits gets the same lower-order contribution.
        block = 1 << (shift * (q - r - 1))
        out[:, 0] += np.repeat(lo.reshape(-1), block)[: rows]

    return _to_lists(_rescale(out, log_base))


# --------------------------------------------------------------- score range
def max_score(mat: Matrix, a: int | None = None) -> float:
    """Highest score the matrix can assign to any sequence."""
    return _extreme_score(mat, a, best=True)


def min_score(mat: Matrix, a: int | None = None) -> float:
    """Lowest score the matrix can assign to any sequence."""
    return _extreme_score(mat, a, best=False)


def _extreme_score(mat: Matrix, a: int | None, best: bool) -> float:
    matrix = _as_array(mat)
    rows, cols = matrix.shape
    if a is None:
        column = matrix.max(axis=0) if best else matrix.min(axis=0)
        # Accumulate left to right rather than with ndarray.sum(), whose pairwise
        # summation would land a few ULP away from the C++ implementation and
        # shift the thresholds that fall back on these bounds.
        total = 0.0
        for value in column:
            total += float(value)
        return total

    # High-order: the score of a column depends on the preceding q-1 symbols, so
    # carry the best total reachable for each possible context.
    q = misc.q_gram_size(rows, a)
    shift = misc.shift(a)
    context_count = 1 << (shift * (q - 1))
    context_mask = context_count - 1

    rows_index = np.arange(rows)
    incoming = rows_index >> shift  # context before this symbol
    outgoing = rows_index & context_mask  # context after it

    pick = np.maximum if best else np.minimum
    fill = -np.inf if best else np.inf

    totals = np.zeros(context_count)
    for i in range(cols):
        candidates = matrix[:, i] + totals[incoming]
        updated = np.full(context_count, fill)
        pick.at(updated, outgoing, candidates)
        totals = updated
    return float(totals.max() if best else totals.min())


def min_delta(mat: Matrix) -> float:
    """The narrowest gap between the best and second-best symbol of any column.

    Symbols tied with the best count as one, so a column whose top two scores
    are equal contributes the gap to the next distinct value below them.
    """
    matrix = _as_array(mat)
    smallest = math.inf
    for column in matrix.T:
        best = second = -math.inf
        for value in column:
            if value > best:
                second = best
                best = value
            elif second < value < best:
                second = value
        smallest = min(smallest, best - second)
    return float(smallest)


# ------------------------------------------------------------------ p-values
def threshold_from_p(
    pssm: Matrix,
    bg: Sequence[float],
    p: float,
    a: int | None = None,
    precision: float | None = None,
) -> float:
    """Lowest score whose false-positive rate under ``bg`` still exceeds ``p``.

    Scores are quantised to multiples of ``1 / precision`` and their
    distribution under ``bg`` computed exactly by dynamic programming, so the
    answer is only as fine-grained as ``precision`` allows. Pass ``a`` for a
    high-order matrix.
    """
    matrix = _as_array(pssm)
    rows = matrix.shape[0]
    if not 0 < p <= 1:
        raise ValueError(f"p-value must be in (0, 1], got {p}")
    bg_arr = np.array(check_bg(bg, rows if a is None else a))
    if precision is None:
        precision = DEFAULT_DP_PRECISION
    precision = float(precision)
    if precision <= 0:
        raise ValueError(f"precision must be positive, got {precision}")

    quantised, min_value, span = _quantise(matrix, precision)
    cols = matrix.shape[1]

    if a is None:
        tail = _score_distribution(quantised, bg_arr, min_value, span)
    else:
        tail = _score_distribution_high_order(
            quantised, bg_arr, min_value, span, a, rows
        )

    return _threshold_from_distribution(
        tail, span, p, precision, cols, min_value, matrix, a
    )


def threshold_from_p_with_precision(
    pssm: Matrix,
    bg: Sequence[float],
    p: float,
    precision: float,
    a: int | None = None,
) -> float:
    """:func:`threshold_from_p` with the MOODS argument order."""
    return threshold_from_p(pssm, bg, p, a=a, precision=precision)


def _quantise(matrix: np.ndarray, precision: float) -> tuple[np.ndarray, int, int]:
    """Scale scores to integers, rounding halves away from zero.

    Returns the integer matrix, the smallest per-column minimum (the offset that
    makes every quantised score non-negative once subtracted), and the width of
    the resulting score range.
    """
    scaled = precision * matrix
    quantised = np.trunc(scaled + np.where(matrix > 0.0, 0.5, -0.5)).astype(np.int64)
    cols = matrix.shape[1]
    min_value = int(quantised.min(axis=0).min())
    span = int(quantised.max(axis=0).sum()) - cols * min_value
    return quantised, min_value, span


def _score_distribution(
    quantised: np.ndarray, bg: np.ndarray, min_value: int, span: int
) -> np.ndarray:
    """Probability of each attainable total score, for a 0-order matrix."""
    rows, cols = quantised.shape
    table = np.zeros(span + 1)
    for j in range(rows):
        table[quantised[j, 0] - min_value] += bg[j]

    for i in range(1, cols):
        nxt = np.zeros(span + 1)
        for j in range(rows):
            s = int(quantised[j, i]) - min_value
            nxt[s:] += bg[j] * table[: span + 1 - s]
        table = nxt
    return table


def _score_distribution_high_order(
    quantised: np.ndarray,
    bg: np.ndarray,
    min_value: int,
    span: int,
    a: int,
    rows: int,
) -> np.ndarray:
    """Probability of each attainable total score, for a high-order matrix.

    The state carried between columns is the trailing ``q - 1`` symbols, since
    those are the context the next column's score depends on.
    """
    cols = quantised.shape[1]
    q = misc.q_gram_size(rows, a)
    shift = misc.shift(a)
    symbol_mask = (1 << shift) - 1
    context_count = 1 << (shift * (q - 1))
    context_mask = context_count - 1

    table = np.zeros((context_count, span + 1))
    for code in range(rows):
        probability = 1.0
        for i in range(q):
            probability *= bg[(code >> (shift * (q - i - 1))) & symbol_mask]
        table[code & context_mask, quantised[code, 0] - min_value] += probability

    for i in range(1, cols):
        nxt = np.zeros((context_count, span + 1))
        for code in range(rows):
            before = (code >> shift) & context_mask
            after = code & context_mask
            symbol = code & symbol_mask
            s = int(quantised[code, i]) - min_value
            nxt[after, s:] += bg[symbol] * table[before, : span + 1 - s]
        table = nxt
    return table.sum(axis=0)


def _threshold_from_distribution(
    tail: np.ndarray,
    span: int,
    p: float,
    precision: float,
    cols: int,
    min_value: int,
    matrix: np.ndarray,
    a: int | None,
) -> float:
    """Walk down the score distribution until the tail probability passes ``p``."""
    cumulative = np.cumsum(tail[::-1])

    # The whole top bucket already exceeds p: no attainable score is rare enough,
    # so sit just below the best score, halfway into the narrowest column gap.
    if cumulative[0] > p:
        return max_score(matrix, a) - min_delta(matrix) / 2

    exceeded = np.flatnonzero(cumulative > p)
    if exceeded.size == 0:
        # Even the entire distribution stays under p; nothing can be excluded.
        return min_score(matrix, a) - 1.0

    r = span - int(exceeded[0])
    return (r + cols * min_value + 1) / precision


# ------------------------------------------------------------------- helpers
def _as_array(mat: Matrix) -> np.ndarray:
    """Validate a matrix and return it as a 2-D float array."""
    if isinstance(mat, np.ndarray) and mat.ndim == 2 and mat.dtype == np.float64:
        if mat.size == 0:
            raise ValueError("matrix has no rows or no columns")
        return mat
    rows = list(mat)
    if not rows:
        raise ValueError("matrix has no rows")
    width = len(rows[0])
    if width == 0:
        raise ValueError("matrix has no columns")
    if any(len(row) != width for row in rows):
        raise ValueError("matrix rows have differing lengths")
    return np.asarray(rows, dtype=np.float64)


def _log(values: np.ndarray) -> np.ndarray:
    """Natural logarithm, elementwise through :func:`math.log`.

    numpy's vectorised ``log`` can land one ULP away from the platform's libm,
    which is the one the C++ implementation calls. Matrices are small enough
    that going through ``math.log`` costs nothing and keeps the two bit-identical.
    Zero maps to ``-inf`` and negatives to ``nan``, as ``np.log`` would.
    """
    flat = np.asarray(values, dtype=np.float64).ravel()
    out = np.fromiter(
        (
            math.log(v) if v > 0.0 else (-math.inf if v == 0.0 else math.nan)
            for v in flat
        ),
        dtype=np.float64,
        count=flat.size,
    )
    return out.reshape(np.shape(values))


def _to_lists(matrix: np.ndarray) -> list[list[float]]:
    return [[float(v) for v in row] for row in matrix]


def _rescale(scores: np.ndarray, log_base: float | None) -> np.ndarray:
    if log_base is None:
        return scores
    if log_base <= 1:
        raise ValueError(f"log_base must be > 1, got {log_base}")
    return scores / math.log(log_base)


def _is_matrix(value) -> bool:
    """True for a sequence of sequences, i.e. a matrix rather than a background."""
    try:
        first = value[0]
    except (TypeError, IndexError, KeyError):
        return False
    return isinstance(first, (list, tuple, np.ndarray)) or (
        hasattr(first, "__len__") and not isinstance(first, (str, bytes))
    )
