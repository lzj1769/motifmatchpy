"""The scanning inner loops, compiled with numba.

Everything here works on flat numpy arrays rather than objects, because that is
what numba can compile to machine code. :mod:`motifmatchpy.scan` packs the
motifs into these arrays once and then reuses them for every sequence.

The algorithm is MOODS': slide a window of ``l`` symbols along the sequence,
pack its contents into an integer, and look that integer up in a table listing
the motifs that could still reach their threshold. Motifs no longer than the
window are decided by the lookup alone; longer ones get scored the rest of the
way by :func:`check_hit_0` or :func:`check_hit_high`, which give up as soon as
the best possible completion falls short.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

try:  # pragma: no cover - exercised by whichever branch the environment takes
    from numba import njit
except ImportError:  # pragma: no cover
    def njit(*args, **kwargs):
        """Fallback so the package imports (slowly) without numba installed."""
        def wrap(func):
            return func

        return wrap(args[0]) if args and callable(args[0]) else wrap


__all__ = ["MotifTables", "process_matches", "ALL_HITS", "MAX_HITS", "COUNTS"]

#: Collect every match.
ALL_HITS = 0
#: Collect matches, but stop recording a motif once it reaches ``max_hits``.
MAX_HITS = 1
#: Count matches only, with the same per-motif cut-off.
COUNTS = 2

_INITIAL_CAPACITY = 1024


class MotifTables(NamedTuple):
    """Every motif packed into flat arrays, ready for :func:`process_matches`.

    ``hit_*`` is the window lookup table in compressed-row form: the motifs
    listed for window code ``c`` are ``hit_index[c] : hit_index[c + 1]``.

    The per-motif matrices live end to end in ``mat_data``, row-major, starting
    at ``mat_index[k]``; the lookahead tables and the high-order prefix/suffix
    tables are laid out the same way.
    """

    # window lookup table
    hit_index: np.ndarray
    hit_motif: np.ndarray
    hit_score: np.ndarray
    hit_full: np.ndarray
    # per motif
    kind: np.ndarray
    size: np.ndarray
    window_pos: np.ndarray
    threshold: np.ndarray
    # score matrices
    mat_index: np.ndarray
    mat_data: np.ndarray
    mat_cols: np.ndarray
    # 0-order lookahead
    la_index: np.ndarray
    la_order: np.ndarray
    la_score: np.ndarray
    # high-order bit layout and reachability tables
    q: np.ndarray
    shift: np.ndarray
    gram_mask: np.ndarray
    ctx_shift: np.ndarray
    ctx_mask: np.ndarray
    ctx_size: np.ndarray
    prefix_index: np.ndarray
    prefix_data: np.ndarray
    suffix_index: np.ndarray
    suffix_data: np.ndarray


@njit(cache=True)
def check_hit_0(
    k, codes, window_match_pos, score, l,
    size, window_pos, threshold, mat_index, mat_data, mat_cols,
    la_index, la_order, la_score,
):
    """Finish scoring a 0-order motif whose window already looked promising.

    Positions outside the window are visited in the lookahead order, most
    discriminating first, and the moment the score plus everything still
    achievable drops below the threshold there is no point continuing.
    """
    m = size[k]
    if m <= l:
        return True, score

    wp = window_pos[k]
    limit = threshold[k]
    origin = window_match_pos - wp
    base = mat_index[k]
    cols = mat_cols[k]
    la_base = la_index[k]

    for i in range(m - l):
        if score + la_score[la_base + i] < limit:
            return False, -np.inf
        column = la_order[la_base + i]
        symbol = np.int64(codes[origin + column])
        score += mat_data[base + symbol * cols + column]
    return score >= limit, score


@njit(cache=True)
def check_hit_high(
    k, codes, window_match_pos, score, l,
    size, window_pos, threshold, mat_index, mat_data, mat_cols,
    q, shift, gram_mask, ctx_shift, ctx_mask, ctx_size,
    prefix_index, prefix_data, suffix_index, suffix_data,
):
    """Finish scoring a high-order motif, working outwards from the window.

    The part before the window is scored backwards and the part after it
    forwards, each abandoned as soon as the best remaining prefix or suffix
    score cannot make up the shortfall.
    """
    m = size[k]
    if m <= l:
        return True, score

    wp = window_pos[k]
    limit = threshold[k]
    origin = window_match_pos - wp
    base = mat_index[k]
    cols = mat_cols[k]
    qk = q[k]
    sh = shift[k]
    gmask = gram_mask[k]
    cshift = ctx_shift[k]
    cmask = ctx_mask[k]
    csize = ctx_size[k]
    has_suffix = wp < m - l

    # Context spanning the window's trailing edge, needed by both halves.
    back_code = np.int64(0)
    if has_suffix:
        for i in range(qk - 1):
            back_code = (back_code << sh) ^ np.int64(codes[origin + wp + l - qk + i + 1])

    if wp > 0:
        forward_limit = limit
        if has_suffix:
            forward_limit -= suffix_data[suffix_index[k] + back_code]

        code = np.int64(0)
        for i in range(qk):
            code = (code << sh) ^ np.int64(codes[origin + wp + i - 1])
        score += mat_data[base + code * cols + (wp - 1)]

        for i in range(1, wp):
            reachable = prefix_data[prefix_index[k] + (wp - i - 1) * csize + (code >> sh)]
            if reachable + score < forward_limit:
                return False, score
            code = (code >> sh) ^ (np.int64(codes[origin + wp - i - 1]) << cshift)
            score += mat_data[base + code * cols + (wp - i - 1)]

    if has_suffix:
        code = back_code
        first = wp + l - qk + 1
        for i in range(first, cols):
            reachable = suffix_data[suffix_index[k] + (i - first) * csize + (code & cmask)]
            if reachable + score < limit:
                return False, score
            code = (gmask & (code << sh)) ^ np.int64(codes[origin + i + qk - 1])
            score += mat_data[base + code * cols + i]

    return score >= limit, score


@njit(cache=True)
def process_matches(codes, bounds, shift, window_size, tables, mode, max_hits):
    """Scan ``codes`` and report every match above threshold.

    ``bounds`` holds the regions to scan as alternating start/end offsets;
    anything outside them contains symbols that are not in the alphabet, and a
    match may not span one.

    Returns ``(motif, pos, score, counts)``. In :data:`COUNTS` mode only
    ``counts`` is filled.
    """
    n_motifs = tables.size.shape[0]
    window_mask = (np.int64(1) << np.int64(shift * window_size)) - 1
    l = window_size

    counts = np.zeros(n_motifs, dtype=np.int64)
    active = np.ones(n_motifs, dtype=np.bool_)

    capacity = _INITIAL_CAPACITY
    out_motif = np.empty(capacity, dtype=np.int64)
    out_pos = np.empty(capacity, dtype=np.int64)
    out_score = np.empty(capacity, dtype=np.float64)
    found = 0

    limited = mode != ALL_HITS
    store = mode != COUNTS

    for region in range(0, bounds.shape[0], 2):
        start = bounds[region]
        end = bounds[region + 1]

        if end - start < l:
            # The region is shorter than the window, so every window that
            # touches it is padded at the tail; only motifs that fit entirely in
            # what is left of the region can match.
            code = np.int64(0)
            for i in range(start, end):
                code = (code << shift) + np.int64(codes[i])
            for _ in range(end - start, l - 1):
                code = (code << shift) & window_mask

            for i in range(start, end):
                code = (code << shift) & window_mask
                for h in range(tables.hit_index[code], tables.hit_index[code + 1]):
                    k = tables.hit_motif[h]
                    if not active[k]:
                        continue
                    if tables.hit_full[h] == 1 and tables.size[k] <= end - i:
                        counts[k] += 1
                        if store:
                            if found == capacity:
                                capacity *= 2
                                out_motif = _grow_int(out_motif, capacity, found)
                                out_pos = _grow_int(out_pos, capacity, found)
                                out_score = _grow_float(out_score, capacity, found)
                            out_motif[found] = k
                            out_pos[found] = i
                            out_score[found] = tables.hit_score[h]
                            found += 1
                        if limited and counts[k] >= max_hits:
                            active[k] = False
            continue

        # Fill the window with the first l - 1 symbols; the loop below adds the
        # last one and then keeps sliding.
        code = np.int64(0)
        for i in range(start, start + l - 1):
            code = (code << shift) + np.int64(codes[i])

        for i in range(start, end - l + 1):
            code = ((code << shift) + np.int64(codes[i + l - 1])) & window_mask
            for h in range(tables.hit_index[code], tables.hit_index[code + 1]):
                k = tables.hit_motif[h]
                if not active[k]:
                    continue

                if tables.hit_full[h] == 1:
                    # The motif fits inside the window, so the lookup settled it.
                    pos = i
                    score = tables.hit_score[h]
                else:
                    wp = tables.window_pos[k]
                    # Skip candidates the motif would overhang the region for.
                    if i - start < wp or i + tables.size[k] - wp > end:
                        continue
                    if tables.kind[k] == 0:
                        matched, score = check_hit_0(
                            k, codes, i, tables.hit_score[h], l,
                            tables.size, tables.window_pos, tables.threshold,
                            tables.mat_index, tables.mat_data, tables.mat_cols,
                            tables.la_index, tables.la_order, tables.la_score,
                        )
                    else:
                        matched, score = check_hit_high(
                            k, codes, i, tables.hit_score[h], l,
                            tables.size, tables.window_pos, tables.threshold,
                            tables.mat_index, tables.mat_data, tables.mat_cols,
                            tables.q, tables.shift, tables.gram_mask,
                            tables.ctx_shift, tables.ctx_mask, tables.ctx_size,
                            tables.prefix_index, tables.prefix_data,
                            tables.suffix_index, tables.suffix_data,
                        )
                    if not matched:
                        continue
                    pos = i - wp

                counts[k] += 1
                if store:
                    if found == capacity:
                        capacity *= 2
                        out_motif = _grow_int(out_motif, capacity, found)
                        out_pos = _grow_int(out_pos, capacity, found)
                        out_score = _grow_float(out_score, capacity, found)
                    out_motif[found] = k
                    out_pos[found] = pos
                    out_score[found] = score
                    found += 1
                if limited and counts[k] >= max_hits:
                    active[k] = False

        # Windows that run off the end of the region are padded, so only motifs
        # short enough to still fit inside the region can match there.
        for i in range(end - l + 1, end):
            code = (code << shift) & window_mask
            for h in range(tables.hit_index[code], tables.hit_index[code + 1]):
                k = tables.hit_motif[h]
                if not active[k]:
                    continue
                if tables.hit_full[h] == 1 and tables.size[k] <= end - i:
                    counts[k] += 1
                    if store:
                        if found == capacity:
                            capacity *= 2
                            out_motif = _grow_int(out_motif, capacity, found)
                            out_pos = _grow_int(out_pos, capacity, found)
                            out_score = _grow_float(out_score, capacity, found)
                        out_motif[found] = k
                        out_pos[found] = i
                        out_score[found] = tables.hit_score[h]
                        found += 1
                    if limited and counts[k] >= max_hits:
                        active[k] = False

    return out_motif[:found], out_pos[:found], out_score[:found], counts


@njit(cache=True)
def _grow_int(values, capacity, used):
    out = np.empty(capacity, dtype=np.int64)
    out[:used] = values[:used]
    return out


@njit(cache=True)
def _grow_float(values, capacity, used):
    out = np.empty(capacity, dtype=np.float64)
    out[:used] = values[:used]
    return out
