"""Agreement with the original MOODS C++ implementation.

Every test here runs the same input through ``motifmatchpy`` and through the
C++ core in ``reference/``, and asserts the answers are *identical* -- bit for
bit for the numeric ones, not merely close. Where the two deliberately differ,
because the C++ has a bug that would be wrong to reproduce, the difference is
recorded in ``reference/core/VENDORED.md`` and patched in that copy, so these
tests still compare like with like.

The whole module is skipped when the reference build is not installed.
"""

from __future__ import annotations

import itertools

import pytest

from motifmatchpy import parsers, scan, tools

moods = pytest.importorskip("moods_reference")

pytestmark = pytest.mark.reference

BACKGROUNDS = [
    [0.25, 0.25, 0.25, 0.25],
    [0.31, 0.19, 0.23, 0.27],
    [0.40, 0.10, 0.10, 0.40],
    [0.70, 0.10, 0.10, 0.10],
]

P_VALUES = [1e-6, 1e-4, 1e-2, 0.3, 1.0]
PRECISIONS = [None, 100.0, 2000.0, 20000.0]


def matches(results):
    """Positions and scores, comparable across the two implementations."""
    return [(m.pos, m.score) for m in results]


def variant_matches(results):
    return [(m.pos, m.score, tuple(m.variants)) for m in results]


def flat(matrix):
    return [float(v) for row in matrix for v in row]


# ------------------------------------------------------------------- tools
@pytest.mark.parametrize("alphabet_size", [1, 2, 4, 5, 20])
def test_flat_bg(alphabet_size):
    assert tools.flat_bg(alphabet_size) == moods.tools.flat_bg(alphabet_size)


@pytest.mark.parametrize("ps", [0.0, 0.01, 1.0, 7.5])
def test_bg_from_sequence(dna, ps):
    assert tools.bg_from_sequence_dna(dna, ps) == moods.tools.bg_from_sequence_dna(
        dna, ps
    )


def test_snp_variants(dna):
    seq = "".join("W" if i % 37 == 0 else c for i, c in enumerate(dna[:5000]))
    mine = [(v.start_pos, v.end_pos, v.modified_seq) for v in tools.snp_variants(seq)]
    theirs = [
        (v.start_pos, v.end_pos, v.modified_seq) for v in moods.tools.snp_variants(seq)
    ]
    assert mine == theirs


def test_snp_variants_for_every_ambiguity_code():
    seq = "ACGTWSMKRYBDHVwsmkrybdhvNn"
    mine = [(v.start_pos, v.end_pos, v.modified_seq) for v in tools.snp_variants(seq)]
    theirs = [
        (v.start_pos, v.end_pos, v.modified_seq) for v in moods.tools.snp_variants(seq)
    ]
    assert mine == theirs


@pytest.mark.parametrize("bg", BACKGROUNDS)
@pytest.mark.parametrize("ps", [0.0, 0.01, 1.0])
@pytest.mark.parametrize("log_base", [None, 2.0, 10.0])
def test_log_odds(pfm_files, bg, ps, log_base):
    for path in pfm_files:
        counts = parsers.pfm(path)
        mine = tools.log_odds(counts, bg, ps, log_base)
        theirs = (
            moods.tools.log_odds(counts, bg, ps)
            if log_base is None
            else moods.tools.log_odds(counts, bg, ps, log_base)
        )
        assert flat(mine) == flat(theirs), path.name


@pytest.mark.parametrize("bg", BACKGROUNDS)
@pytest.mark.parametrize("ps", [0.0, 0.01, 1.0])
def test_log_odds_high_order(adm_file, bg, ps):
    first = parsers.adm_1o_terms(adm_file)
    low = [[row[0] for row in parsers.adm_0o_terms(adm_file)]]
    mine = tools.log_odds_high_order(first, low, bg, ps, 4)
    theirs = moods.tools.log_odds(first, low, bg, ps, 4)
    assert flat(mine) == flat(theirs)


def test_reverse_complement(pfm_files, adm_file, flat4):
    for path in pfm_files:
        matrix = parsers.pfm_to_log_odds(path, flat4, 0.01)
        assert flat(tools.reverse_complement(matrix)) == flat(
            moods.tools.reverse_complement(matrix)
        )
        assert flat(tools.reverse_complement(matrix, 4)) == flat(
            moods.tools.reverse_complement(matrix, 4)
        )
    high = parsers.adm_to_log_odds(adm_file, flat4, 0.01)
    assert flat(tools.reverse_complement(high, 4)) == flat(
        moods.tools.reverse_complement(high, 4)
    )


@pytest.mark.parametrize("bg", BACKGROUNDS)
def test_score_range(pfm_files, adm_file, bg):
    for path in pfm_files:
        matrix = parsers.pfm_to_log_odds(path, bg, 0.01)
        assert tools.max_score(matrix) == moods.tools.max_score(matrix)
        assert tools.min_score(matrix) == moods.tools.min_score(matrix)
        assert tools.max_score(matrix, 4) == moods.tools.max_score(matrix, 4)
        assert tools.min_score(matrix, 4) == moods.tools.min_score(matrix, 4)
        assert tools.min_delta(matrix) == moods.tools.min_delta(matrix)

    high = parsers.adm_to_log_odds(adm_file, bg, 0.01)
    assert tools.max_score(high, 4) == moods.tools.max_score(high, 4)
    assert tools.min_score(high, 4) == moods.tools.min_score(high, 4)


@pytest.mark.parametrize("bg", BACKGROUNDS)
@pytest.mark.parametrize("p", P_VALUES)
@pytest.mark.parametrize("precision", PRECISIONS)
def test_threshold_from_p(pfm_files, bg, p, precision):
    for path in pfm_files:
        matrix = parsers.pfm_to_log_odds(path, bg, 0.01)
        for a in (None, 4):
            mine = tools.threshold_from_p(matrix, bg, p, a=a, precision=precision)
            if precision is None:
                theirs = (
                    moods.tools.threshold_from_p(matrix, bg, p)
                    if a is None
                    else moods.tools.threshold_from_p(matrix, bg, p, 4)
                )
            else:
                theirs = (
                    moods.tools.threshold_from_p_with_precision(matrix, bg, p, precision)
                    if a is None
                    else moods.tools.threshold_from_p_with_precision(
                        matrix, bg, p, precision, 4
                    )
                )
            assert mine == theirs, f"{path.name} a={a}"


@pytest.mark.parametrize("bg", BACKGROUNDS)
@pytest.mark.parametrize("p", P_VALUES)
@pytest.mark.parametrize("precision", PRECISIONS)
def test_threshold_from_p_high_order(adm_file, bg, p, precision):
    matrix = parsers.adm_to_log_odds(adm_file, bg, 0.01)
    mine = tools.threshold_from_p(matrix, bg, p, a=4, precision=precision)
    theirs = (
        moods.tools.threshold_from_p(matrix, bg, p, 4)
        if precision is None
        else moods.tools.threshold_from_p_with_precision(matrix, bg, p, precision, 4)
    )
    assert mine == theirs


# ----------------------------------------------------------------- parsers
@pytest.mark.parametrize("bg", BACKGROUNDS)
@pytest.mark.parametrize("log_base", [None, 2.0])
def test_parsers(pfm_files, adm_file, bg, log_base):
    base = -1 if log_base is None else log_base
    for path in pfm_files:
        assert flat(parsers.pfm(path)) == flat(moods.parsers.pfm(str(path)))
        assert flat(parsers.pfm_to_log_odds(path, bg, 0.01, log_base)) == flat(
            moods.parsers.pfm_to_log_odds(str(path), bg, 0.01, base)
        )
    assert flat(parsers.adm_1o_terms(adm_file)) == flat(
        moods.parsers.adm_1o_terms(str(adm_file), 4)
    )
    assert flat(parsers.adm_0o_terms(adm_file)) == flat(
        moods.parsers.adm_0o_terms(str(adm_file), 4)
    )
    assert flat(parsers.adm_to_log_odds(adm_file, bg, 0.01, 4, log_base)) == flat(
        moods.parsers.adm_to_log_odds(str(adm_file), bg, 0.01, 4, base)
    )


# -------------------------------------------------------------------- misc
def test_misc_primitives():
    from motifmatchpy import misc
    from motifmatchpy._alphabet import build_map

    for a in (2, 4, 5, 8, 16):
        assert misc.shift(a) == moods.misc.shift(a)
        assert misc.mask(a) == moods.misc.mask(a)
    for rows in (4, 16, 64):
        assert misc.q_gram_size(rows, 4) == moods.misc.q_gram_size(rows, 4)
    for code in range(64):
        assert misc.rc_tuple(code, 4, 3) == moods.misc.rc_tuple(code, 4, 3)

    table = [int(x) for x in build_map()]
    for seq in ("", "ACGT", "NNACGTNN", "acgtNNacgt", "N" * 10):
        assert misc.preprocess_seq(seq, 4, table) == list(
            moods.misc.preprocess_seq(seq, 4, table)
        )


# ----------------------------------------------------------------- scanning
@pytest.mark.parametrize("bg", BACKGROUNDS[:2])
@pytest.mark.parametrize("p", [1e-5, 1e-4, 1e-3])
@pytest.mark.parametrize("window_size", [4, 7, 8])
def test_scan_dna(dna, pfm_files, bg, p, window_size):
    matrices = [parsers.pfm_to_log_odds(path, bg, 0.01) for path in pfm_files]
    thresholds = [tools.threshold_from_p(m, bg, p, a=4) for m in matrices]
    mine = scan.scan_dna(dna, matrices, bg, thresholds, window_size)
    theirs = moods.scan.scan_dna(dna, matrices, bg, thresholds, window_size)
    assert [matches(r) for r in mine] == [matches(r) for r in theirs]


@pytest.mark.parametrize("bg", BACKGROUNDS[:2])
@pytest.mark.parametrize("window_size", [3, 5, 7])
def test_scan_dna_high_order(dna, adm_file, bg, window_size):
    matrix = parsers.adm_to_log_odds(adm_file, bg, 0.01)
    threshold = tools.threshold_from_p(matrix, bg, 1e-4, a=4)
    mine = scan.scan_dna(dna, [matrix], bg, [threshold], window_size)
    theirs = moods.scan.scan_dna(dna, [matrix], bg, [threshold], window_size)
    assert matches(mine[0]) == matches(theirs[0])


def test_scan_dna_mixed_orders_and_strands(dna, pfm_files, adm_file, flat4):
    matrices = [parsers.pfm_to_log_odds(p, flat4, 0.01) for p in pfm_files]
    matrices.append(parsers.adm_to_log_odds(adm_file, flat4, 0.01))
    matrices += [tools.reverse_complement(m, 4) for m in matrices]
    thresholds = [tools.threshold_from_p(m, flat4, 1e-4, a=4) for m in matrices]
    mine = scan.scan_dna(dna, matrices, flat4, thresholds)
    theirs = moods.scan.scan_dna(dna, matrices, flat4, thresholds)
    assert [matches(r) for r in mine] == [matches(r) for r in theirs]


def test_custom_alphabet_scan(rng):
    seq = "".join(rng.choice("ACGTM") for _ in range(20000))
    matrix = [
        [0.5, -1.0, 0.3],
        [0.1, 0.2, -0.4],
        [-0.3, 0.4, 0.8],
        [0.7, -0.2, 0.1],
        [1.0, 1.0, -0.5],
    ]
    alphabet = ["aA", "cC", "gG", "tT", "mM"]
    bg = tools.flat_bg(5)
    for window_size in (3, 5, 6):
        mine = scan.scan(seq, [matrix], bg, [1.0], window_size, alphabet)
        theirs = moods.scan.scan(seq, [matrix], bg, [1.0], window_size, alphabet)
        assert matches(mine[0]) == matches(theirs[0])


@pytest.mark.parametrize("max_hits", [0, 1, 5, 50, 1000])
def test_max_hits_and_counts(dna, pfm_files, flat4, max_hits):
    matrices = [parsers.pfm_to_log_odds(p, flat4, 0.01) for p in pfm_files]
    thresholds = [tools.threshold_from_p(m, flat4, 1e-2, a=4) for m in matrices]

    mine = scan.Scanner(7)
    mine.set_motifs(matrices, flat4, thresholds)
    theirs = moods.scan.Scanner(7)
    theirs.set_motifs(matrices, flat4, thresholds)

    assert [matches(r) for r in mine.scan_max_hits(dna, max_hits)] == [
        matches(r) for r in theirs.scan_max_hits(dna, max_hits)
    ]
    assert mine.counts_max_hits(dna, max_hits) == list(
        theirs.counts_max_hits(dna, max_hits)
    )


@pytest.mark.parametrize("target", [10, 100, 1000])
@pytest.mark.parametrize("params", [(10, 2, 10), (0, 3, 10), (3, 2, 5)])
def test_best_hits(dna, pfm_files, flat4, target, params):
    iterations, mult, limit = params
    matrices = [parsers.pfm_to_log_odds(p, flat4, 0.01) for p in pfm_files]
    mine = scan.scan_best_hits_dna(dna, matrices, target, iterations, mult, limit)
    theirs = moods.scan.scan_best_hits_dna(
        dna, matrices, target, iterations, mult, limit
    )
    assert [matches(r) for r in mine] == [matches(r) for r in theirs]


def test_naive_scan(dna, pfm_file, adm_file, flat4):
    sample = dna[:50000]
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-3)
    assert matches(scan.naive_scan_dna(sample, matrix, threshold)) == matches(
        moods.scan.naive_scan_dna(sample, matrix, threshold)
    )

    high = parsers.adm_to_log_odds(adm_file, flat4, 0.01)
    high_threshold = tools.threshold_from_p(high, flat4, 1e-3, a=4)
    assert matches(scan.naive_scan_dna(sample, high, high_threshold, 4)) == matches(
        moods.scan.naive_scan_dna(sample, high, high_threshold, 4)
    )


# ----------------------------------------------------------------- variants
def test_variant_matches_for_snps(rng, polya, flat4):
    mine = scan.Scanner(7)
    mine.set_motifs([polya], flat4, [5.0])
    theirs = moods.scan.Scanner(7)
    theirs.set_motifs([polya], flat4, [5.0])

    for _ in range(200):
        seq = "".join(
            rng.choice("ACGTACGTWSMKRYBDHV") for _ in range(rng.randint(6, 40))
        )
        my_variants = tools.snp_variants(seq)
        their_variants = moods.tools.snp_variants(seq)
        for depth in (0, 1, 2):
            assert variant_matches(
                mine.variant_matches(seq, my_variants, depth)[0]
            ) == variant_matches(theirs.variant_matches(seq, their_variants, depth)[0])


def test_variant_matches_for_indels(rng, polya, flat4):
    mine = scan.Scanner(7)
    mine.set_motifs([polya], flat4, [5.0])
    theirs = moods.scan.Scanner(7)
    theirs.set_motifs([polya], flat4, [5.0])

    for _ in range(200):
        length = rng.randint(10, 40)
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        spec = []
        for _ in range(rng.randint(1, 4)):
            start = rng.randint(0, length - 1)
            end = min(length, start + rng.randint(0, 3))
            replacement = "".join(
                rng.choice("ACGT") for _ in range(rng.randint(0, 3))
            )
            spec.append((start, end, replacement))
        my_variants = [tools.Variant(*v) for v in spec]
        their_variants = [moods.Variant(*v) for v in spec]
        for depth in (0, 1, 3):
            assert variant_matches(
                mine.variant_matches(seq, my_variants, depth)[0]
            ) == variant_matches(
                theirs.variant_matches(seq, their_variants, depth)[0]
            ), f"{seq!r} {spec} depth={depth}"


# ----------------------------------------------------- exhaustive small cases
@pytest.mark.slow
def test_exhaustive_short_sequences(flat4):
    """Every ACGN sequence up to length 8, at several window sizes.

    Short sequences and ``N`` runs are where the window machinery falls back on
    its edge cases, and where the C++ original turned out to have two separate
    off-by-one faults, so enumerating them is worth the seconds it costs.
    """
    matrix = [
        [1.0, 0.2, 0.3, 0.1],
        [0.0, 0.0, 0.0, 0.0],
        [0.5, 0.1, 0.0, 0.7],
        [0.1, 0.8, 0.2, 0.0],
    ]
    for length in range(1, 9):
        for letters in itertools.product("ACGN", repeat=length):
            seq = "".join(letters)
            for window_size in (2, 5, 7):
                mine = scan.scan_dna(seq, [matrix], flat4, [1.5], window_size)
                theirs = moods.scan.scan_dna(seq, [matrix], flat4, [1.5], window_size)
                assert matches(mine[0]) == matches(theirs[0]), (
                    f"{seq!r} window={window_size}"
                )


@pytest.mark.slow
def test_random_matrices_against_random_sequences(rng, flat4):
    """Random matrices explore shapes the curated JASPAR set does not."""
    for _ in range(40):
        width = rng.randint(1, 12)
        matrix = [[rng.uniform(-4, 4) for _ in range(width)] for _ in range(4)]
        seq = "".join(
            rng.choice("ACGTACGTACGTN") for _ in range(rng.randint(1, 500))
        )
        threshold = tools.threshold_from_p(matrix, flat4, rng.choice([1e-3, 0.05, 0.5]))
        for window_size in (2, 4, 7):
            mine = scan.scan_dna(seq, [matrix], flat4, [threshold], window_size)
            theirs = moods.scan.scan_dna(seq, [matrix], flat4, [threshold], window_size)
            assert matches(mine[0]) == matches(theirs[0])


@pytest.mark.slow
def test_random_high_order_matrices(rng, flat4):
    for _ in range(20):
        width = rng.randint(2, 10)
        matrix = [[rng.uniform(-3, 3) for _ in range(width)] for _ in range(16)]
        seq = "".join(rng.choice("ACGTACGTN") for _ in range(rng.randint(2, 400)))
        threshold = tools.threshold_from_p(matrix, flat4, 0.01, a=4)
        for window_size in (2, 4, 7):
            mine = scan.scan_dna(seq, [matrix], flat4, [threshold], window_size)
            theirs = moods.scan.scan_dna(seq, [matrix], flat4, [threshold], window_size)
            assert matches(mine[0]) == matches(theirs[0]), (
                f"width={width} window={window_size}"
            )
