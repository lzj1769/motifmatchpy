"""Per-motif preprocessing for the scanner.

MOODS does not score every position against every matrix. It reads a fixed-width
*window* of the sequence, looks the window's contents up in a table that says
which motifs could still reach their threshold, and only then scores the
candidates in full -- abandoning each as soon as the best possible completion
falls short. These classes build the tables that make that work:

* where to place the window within the motif (:func:`window_position`), chosen
  so the window carries as much of the score as possible;
* for a 0-order motif, the order to examine the remaining positions in and the
  best score still reachable after each of them (the *lookahead*);
* for a high-order motif, the best prefix and suffix scores reachable from each
  context, which serve the same purpose.

The arithmetic mirrors ``core/motif_0.cpp`` and ``core/motif_h.cpp`` operation
for operation, including where upstream's own accumulation order is idiosyncratic,
so that scores come out bit-identical rather than merely close.
"""

from __future__ import annotations

import math

import numpy as np

from . import misc

__all__ = ["Motif0", "MotifH", "build_motif"]


def build_motif(matrix, bg, window_size: int, threshold: float, a: int):
    """A :class:`Motif0` or :class:`MotifH`, whichever the matrix calls for."""
    rows = len(matrix)
    return (
        Motif0(matrix, bg, window_size, threshold)
        if rows == a
        else MotifH(matrix, bg, window_size, threshold, a)
    )


class Motif0:
    """A plain position weight matrix: one row per symbol."""

    kind = 0

    def __init__(self, matrix, bg, window_size: int, threshold: float) -> None:
        self.mat = np.asarray(matrix, dtype=np.float64)
        self.a, self.m = self.mat.shape
        self.l = int(window_size)
        self.T = float(threshold)

        differences = expected_differences(self.mat, np.asarray(bg, dtype=np.float64))
        self.wp = window_position(differences, self.l, self.m)
        self.lookahead_order = _lookahead_order(differences, self.l, self.wp, self.m)
        self.lookahead_scores = _lookahead_scores(self.mat, self.lookahead_order, self.l)

    @property
    def size(self) -> int:
        """Number of sequence positions the motif spans."""
        return self.m

    def window_scores(self, shift: int) -> tuple[np.ndarray, np.ndarray]:
        """Score and candidacy of every possible window content.

        Returns ``(is_candidate, score)`` indexed by the packed window code.
        ``score`` is the score of the part of the motif that lies inside the
        window; ``is_candidate`` says whether the best possible completion of
        that window could still reach the threshold.
        """
        codes = np.arange(1 << (shift * self.l), dtype=np.int64)
        symbol_mask = misc.mask(self.a)

        inside = min(self.m, self.l)
        offset = 0 if self.l >= self.m else self.wp

        score = np.zeros(codes.size, dtype=np.float64)
        valid = np.ones(codes.size, dtype=bool)
        for i in range(inside):
            symbols = (codes >> (shift * (self.l - i - 1))) & symbol_mask
            # A window code can name a symbol outside the alphabet whenever the
            # alphabet size is not a power of two; such a window matches nothing.
            valid &= symbols < self.a
            score += self.mat[np.minimum(symbols, self.a - 1), offset + i]

        if self.l >= self.m:
            candidate = score >= self.T
        else:
            candidate = score + self.lookahead_scores[0] >= self.T
        return candidate & valid, np.where(valid, score, -np.inf)


class MotifH:
    """A high-order matrix: one row per q-gram of symbols."""

    kind = 1

    def __init__(
        self, matrix, bg, window_size: int, threshold: float, alphabet_size: int
    ) -> None:
        self.mat = np.asarray(matrix, dtype=np.float64)
        self.rows, self.cols = self.mat.shape
        self.a = int(alphabet_size)
        self.l = int(window_size)
        self.T = float(threshold)

        self.q = misc.q_gram_size(self.rows, self.a)
        self.m = self.cols + self.q - 1

        if self.l < self.q:
            raise ValueError(
                f"window size {self.l} is smaller than the motif's q-gram "
                f"length {self.q}"
            )

        self.SHIFT = misc.shift(self.a)
        self.MASK = (1 << (self.SHIFT * self.q)) - 1
        self.Q_SHIFT = self.SHIFT * (self.q - 1)
        self.Q_CODE_SIZE = 1 << self.Q_SHIFT
        self.Q_MASK = self.Q_CODE_SIZE - 1

        expected = self._expected_scores(np.asarray(bg, dtype=np.float64))
        self.wp = self._window_position(expected)
        self.P = self._max_scores_forward(0, self.wp)
        self.S = self._max_scores_backward(self.wp + self.l - self.q + 1, self.cols)

    @property
    def size(self) -> int:
        """Number of sequence positions the motif spans."""
        return self.m

    # -- preprocessing --------------------------------------------------
    def _expected_scores(self, bg: np.ndarray) -> np.ndarray:
        """Mean score of each column under the background."""
        symbol_mask = (1 << self.SHIFT) - 1
        weights = np.ones(self.rows, dtype=np.float64)
        for k in range(self.q):
            symbols = (np.arange(self.rows) >> (self.SHIFT * (self.q - 1 - k))) & symbol_mask
            weights *= bg[symbols]

        expected = np.zeros(self.cols, dtype=np.float64)
        for i in range(self.cols):
            total = 0.0
            for j in range(self.rows):
                total += weights[j] * self.mat[j, i]
            expected[i] = total
        return expected

    def _window_position(self, expected: np.ndarray) -> int:
        """Where to place the window, by largest score-above-expectation.

        Reproduces upstream exactly, including that its running expectation is
        assigned rather than accumulated over the initial window -- so it holds
        one column's mean where the sum was clearly intended. The choice only
        steers the search, never which positions match, but reproducing it keeps
        the accumulation order -- and so the last bits of every score -- the same.
        """
        if self.l >= self.m:
            return 0

        window_scores = np.empty(self.m - self.l + 1, dtype=np.float64)
        for start in range(self.m - self.l + 1):
            reachable = self._max_scores_forward(start, start + self.l - self.q + 1)
            window_scores[start] = reachable[-1].max()

        current = 0.0
        for i in range(self.l - self.q + 1):
            current = expected[i]

        best_loss = window_scores[0] - current
        position = 0
        for i in range(1, self.m - self.l + 1):
            current -= expected[i]
            current += expected[i + self.l - self.q]
            if window_scores[i] - current > best_loss:
                position = i
                best_loss = window_scores[i] - current
        return position

    def _max_scores_forward(self, start: int, end: int) -> np.ndarray:
        """Best score reachable over ``[start, end)`` for each trailing context."""
        width = max(0, end - start)
        scores = np.zeros((width, self.Q_CODE_SIZE), dtype=np.float64)
        if width == 0:
            return scores

        for j in range(self.rows):
            after = j & self.Q_MASK
            scores[0, after] = max(self.mat[j, start], scores[0, after])
        for i in range(1, width):
            for j in range(self.rows):
                after = j & self.Q_MASK
                candidate = self.mat[j, i + start] + scores[i - 1, j >> self.SHIFT]
                scores[i, after] = max(candidate, scores[i, after])
        return scores

    def _max_scores_backward(self, start: int, end: int) -> np.ndarray:
        """Best score reachable over ``[start, end)`` for each leading context."""
        width = max(0, end - start)
        scores = np.zeros((width, self.Q_CODE_SIZE), dtype=np.float64)
        if width == 0:
            return scores

        for j in range(self.rows):
            before = j >> self.SHIFT
            scores[width - 1, before] = max(self.mat[j, end - 1], scores[width - 1, before])
        for i in range(1, width):
            for j in range(self.rows):
                before = j >> self.SHIFT
                candidate = self.mat[j, end - i - 1] + scores[width - i, j & self.Q_MASK]
                scores[width - i - 1, before] = max(
                    candidate, scores[width - i - 1, before]
                )
        return scores

    # -- scanning support -----------------------------------------------
    def window_scores(self, shift: int) -> tuple[np.ndarray, np.ndarray]:
        """Score and candidacy of every possible window content."""
        size = 1 << (shift * self.l)
        candidate = np.zeros(size, dtype=bool)
        score = np.empty(size, dtype=np.float64)

        if self.l >= self.m:
            for code in range(size):
                total = 0.0
                for i in range(self.cols):
                    c = self.MASK & (code >> (self.SHIFT * (self.l - self.q - i)))
                    total += self.mat[c, i]
                score[code] = total
                candidate[code] = total >= self.T
            return candidate, score

        suffix_mask = (1 << (self.SHIFT * (self.q - 1))) - 1
        for code in range(size):
            total = 0.0
            for i in range(self.l - self.q + 1):
                c = self.MASK & (code >> (self.SHIFT * (self.l - self.q - i)))
                total += self.mat[c, self.wp + i]

            potential = total
            if self.wp > 0:
                potential += self.P[-1][code >> (self.SHIFT * (self.l - self.q + 1))]
            if self.wp < self.m - self.l:
                potential += self.S[0][code & suffix_mask]
            score[code] = total
            candidate[code] = potential >= self.T
        return candidate, score


def expected_differences(mat: np.ndarray, bg: np.ndarray) -> np.ndarray:
    """Per column, how far the best symbol beats the background average.

    Columns that discriminate strongly are the ones worth putting inside the
    window, and worth examining first outside it.
    """
    a, m = mat.shape
    out = np.empty(m, dtype=np.float64)
    for i in range(m):
        value = mat[:, i].max()
        for j in range(a):
            value -= bg[j] * mat[j, i]
        out[i] = value
    return out


def window_position(differences: np.ndarray, l: int, m: int) -> int:
    """The window offset whose columns discriminate most."""
    if l >= m:
        return 0

    current = 0.0
    for i in range(l):
        current += differences[i]

    best = current
    position = 0
    for i in range(m - l):
        current -= differences[i]
        current += differences[i + l]
        if current > best:
            best = current
            position = i + 1
    return position


def _lookahead_order(
    differences: np.ndarray, l: int, window_pos: int, m: int
) -> np.ndarray:
    """The positions outside the window, most discriminating first.

    Examining them in this order lets the scan give up as early as possible.
    """
    if l >= m:
        return np.empty(0, dtype=np.int64)
    outside = list(range(window_pos)) + list(range(window_pos + l, m))
    outside.sort(key=lambda i: -differences[i])
    return np.array(outside, dtype=np.int64)


def _lookahead_scores(mat: np.ndarray, order: np.ndarray, l: int) -> np.ndarray:
    """Best score still reachable from each step of the lookahead onwards.

    ``scores[i]`` is the most the positions ``order[i:]`` could contribute, so a
    partial score below ``T - scores[i]`` can never recover and the scan stops.
    """
    if order.size == 0:
        return np.empty(0, dtype=np.float64)

    scores = np.zeros(order.size, dtype=np.float64)
    total = 0.0
    for i in range(order.size - 1, -1, -1):
        best = -math.inf
        for value in mat[:, order[i]]:
            best = max(best, value)
        total += best
        scores[i] = total
    return scores
