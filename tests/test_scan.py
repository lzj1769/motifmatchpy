"""Scanning, checked position by position against a brute-force scan.

``reference.brute_force_scan`` shares no code with the scanner, so agreement
between the two is real evidence that the window/lookahead machinery finds
exactly the right occurrences and scores them correctly.
"""

from __future__ import annotations

import itertools

import pytest

import reference
from motifmatchpy import parsers, scan, tools


def as_pairs(matches, ndigits=9):
    return sorted((m.pos, round(m.score, ndigits)) for m in matches)


@pytest.mark.parametrize("seq, expected", [
    ("ACGT", [(0, 1.0), (1, 0.0), (2, 0.0)]),
    ("", []), ("A", []), ("NNNN", []),
])
@pytest.mark.parametrize("iterations", [0, 1, 10])
def test_best_hits_target_exceeds_available_positions(seq, expected, iterations):
    matrix = [[1.0, 2.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
    found = scan.scan_best_hits_dna(seq, [matrix], 100, iterations=iterations)
    assert [as_pairs(group) for group in found] == [expected]


def rounded(pairs, ndigits=9):
    return sorted((pos, round(score, ndigits)) for pos, score in pairs)


# ------------------------------------------------------- against brute force
def test_scan_dna_finds_exactly_the_brute_force_hits(sample, pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-3)
    found = scan.scan_dna(sample, [matrix], flat4, [threshold])[0]
    expected = reference.brute_force_scan(sample, matrix, threshold)
    assert expected, "the test is vacuous if nothing matches"
    assert as_pairs(found) == rounded(expected)


def test_every_example_matrix_agrees_with_brute_force(dna, pfm_files, flat4):
    sample = dna[:4000]
    matrices = [parsers.pfm_to_log_odds(p, flat4, 0.01) for p in pfm_files]
    thresholds = [tools.threshold_from_p(m, flat4, 1e-2) for m in matrices]
    results = scan.scan_dna(sample, matrices, flat4, thresholds)
    assert len(results) == len(matrices)
    for path, matrix, threshold, found in zip(
        pfm_files, matrices, thresholds, results, strict=True
    ):
        expected = reference.brute_force_scan(sample, matrix, threshold)
        assert as_pairs(found) == rounded(expected), f"mismatch for {path.name}"


def test_high_order_scan_agrees_with_brute_force(sample, adm_file, flat4):
    matrix = parsers.adm_to_log_odds(adm_file, flat4, 0.01)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-3, a=4)
    found = scan.scan_dna(sample, [matrix], flat4, [threshold])[0]
    expected = reference.brute_force_scan(sample, matrix, threshold)
    assert expected
    assert as_pairs(found) == rounded(expected)


def test_mixing_orders_in_one_scan(sample, pfm_file, adm_file, flat4):
    plain = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    high = parsers.adm_to_log_odds(adm_file, flat4, 0.01)
    thresholds = [
        tools.threshold_from_p(plain, flat4, 1e-3, a=4),
        tools.threshold_from_p(high, flat4, 1e-3, a=4),
    ]
    together = scan.scan_dna(sample, [plain, high], flat4, thresholds)
    separately = [
        scan.scan_dna(sample, [plain], flat4, thresholds[:1])[0],
        scan.scan_dna(sample, [high], flat4, thresholds[1:])[0],
    ]
    assert [as_pairs(r) for r in together] == [as_pairs(r) for r in separately]


def test_reverse_complement_hits_match_the_reverse_strand(sample, pfm_file, flat4):
    """A hit for the reverse-complemented matrix at position p must score the
    same as the forward matrix against the reverse complement of that window."""
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    rc = tools.reverse_complement(matrix, 4)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-3)
    hits = scan.scan_dna(sample, [rc], flat4, [threshold])[0]
    assert hits

    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    for hit in hits[:50]:
        window = sample[hit.pos : hit.pos + len(matrix[0])]
        revcomp = window.translate(complement)[::-1]
        assert reference.score_at(revcomp, 0, matrix) == pytest.approx(hit.score)


# ----------------------------------------------------- invalid sequence data
def test_hits_never_span_an_unknown_base(flat4):
    matrix = [[1.0] * 4, [0.0] * 4, [0.0] * 4, [0.0] * 4]  # matches AAAA
    seq = "AAAA" + "N" + "AAAA" + "X" + "AAAA"
    #      0123    4     5678    9     0123
    found = scan.scan_dna(seq, [matrix], flat4, [4.0])[0]
    assert [m.pos for m in found] == [0, 5, 10]


def test_a_sequence_of_only_unknown_bases_yields_nothing(flat4):
    matrix = [[1.0] * 4, [0.0] * 4, [0.0] * 4, [0.0] * 4]
    assert scan.scan_dna("N" * 100, [matrix], flat4, [4.0])[0] == []


def test_an_empty_sequence_yields_nothing(flat4, polya):
    assert scan.scan_dna("", [polya], flat4, [5.0])[0] == []


def test_a_sequence_shorter_than_the_motif_yields_nothing(flat4):
    matrix = [[1.0] * 8, [0.0] * 8, [0.0] * 8, [0.0] * 8]
    assert scan.scan_dna("AAAA", [matrix], flat4, [8.0])[0] == []


def test_a_motif_matching_the_whole_sequence_is_found(flat4):
    matrix = [[1.0] * 4, [0.0] * 4, [0.0] * 4, [0.0] * 4]
    assert [m.pos for m in scan.scan_dna("AAAA", [matrix], flat4, [4.0])[0]] == [0]


def test_lowercase_sequence_scores_the_same(sample, pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-3)
    upper = scan.scan_dna(sample.upper(), [matrix], flat4, [threshold])[0]
    lower = scan.scan_dna(sample.lower(), [matrix], flat4, [threshold])[0]
    assert as_pairs(upper) == as_pairs(lower)


def test_bytes_and_str_give_the_same_hits(sample, pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-3)
    assert as_pairs(scan.scan_dna(sample.encode(), [matrix], flat4, [threshold])[0]) == (
        as_pairs(scan.scan_dna(sample, [matrix], flat4, [threshold])[0])
    )


# ------------------------------------- boundaries: matches at the very edges
def test_hit_ending_at_the_end_of_the_sequence_is_found(flat4):
    """A motif shorter than the window, finishing exactly at the last base.

    The C++ scanner drops these: its end-of-region guard tests
    ``size < end - i`` where the motif fits when ``size <= end - i``.
    """
    matrix = [[1.0] * 4, [0.0] * 4, [0.0] * 4, [0.0] * 4]
    seq = "CCCCCCCCCCAAAA"
    assert [m.pos for m in scan.scan_dna(seq, [matrix], flat4, [4.0], 7)[0]] == [10]


def test_hit_ending_just_before_an_n_is_found(flat4):
    """The same boundary, but at the edge of an ``N`` run rather than the
    sequence end -- which is where it bites on real genomic data."""
    matrix = [[1.0] * 4, [0.0] * 4, [0.0] * 4, [0.0] * 4]
    seq = "CCCCCCAAAANCCCCCC"
    assert [m.pos for m in scan.scan_dna(seq, [matrix], flat4, [4.0], 7)[0]] == [6]


def test_hit_not_at_the_start_of_a_short_region_is_found(flat4):
    """A region shorter than the window, with the match away from its start.

    The C++ scanner loses these: it tests the lookup table with the previous
    position's window code and only advances the window inside that test.
    """
    matrix = [[1.0] * 3, [0.0] * 3, [0.0] * 3, [0.0] * 3]
    assert [m.pos for m in scan.scan_dna("CAAA", [matrix], flat4, [3.0], 7)[0]] == [1]
    assert [m.pos for m in scan.scan_dna("CCAAA", [matrix], flat4, [3.0], 7)[0]] == [2]


def test_high_order_hit_ending_at_the_sequence_end_is_found(flat4):
    matrix = [[0.0] * 3 for _ in range(16)]
    matrix[0][0] = matrix[0][1] = matrix[0][2] = 1.0  # rewards AA at every step
    assert [m.pos for m in scan.scan_dna("CCCCAAAA", [matrix], flat4, [3.0], 7)[0]] == [4]


@pytest.mark.parametrize("window_size", [2, 3, 4, 5, 6, 7, 8])
def test_boundary_hits_are_found_at_every_window_size(flat4, window_size):
    matrix = [[1.0] * 4, [0.0] * 4, [0.0] * 4, [0.0] * 4]
    found = scan.scan_dna("CCCCCCCCCCAAAA", [matrix], flat4, [4.0], window_size)[0]
    assert [m.pos for m in found] == [10]


@pytest.mark.slow
def test_exhaustive_agreement_with_brute_force(flat4):
    """Every ACGN sequence up to length 7, against a scan that shares no code.

    Short sequences and ``N`` runs are where the window machinery has to fall
    back on its edge cases, so enumerating them all is the surest check there is.

    Scores landing exactly on the threshold are excluded: the two implementations
    add the columns up in a different order, so such a hit can fall on either
    side of the cutoff. That is a property of MOODS too -- the C++ scanner agrees
    with this one, not with the brute-force scan, on those cases.
    """
    matrices = [
        [[1.0] * 4, [0.0] * 4, [0.0] * 4, [0.0] * 4],
        [[1.0, 0.5], [0.0, 0.0], [0.2, 0.0], [0.0, 1.0]],
        [
            [1.0, 0.2, 0.3, 0.1, 0.9, 0.4],
            [0.0] * 6,
            [0.5, 0.1, 0.0, 0.7, 0.0, 0.2],
            [0.1, 0.8, 0.2, 0.0, 0.3, 1.0],
        ],
    ]
    for matrix in matrices:
        for length in range(1, 8):
            for letters in itertools.product("ACGN", repeat=length):
                seq = "".join(letters)
                for window_size in (2, 5, 7):
                    found = scan.scan_dna(seq, [matrix], flat4, [1.5], window_size)[0]
                    expected = reference.brute_force_scan(seq, matrix, 1.5)
                    reference.assert_agrees_ignoring_ties(found, expected, 1.5)


# -------------------------------------------------------------- window sizes
@pytest.mark.parametrize("window_size", [1, 2, 4, 7, 8])
def test_window_size_does_not_change_the_answer(sample, pfm_file, flat4, window_size):
    """The window is a search optimisation, so every setting must agree."""
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-3)
    baseline = as_pairs(scan.scan_dna(sample, [matrix], flat4, [threshold], 7)[0])
    assert as_pairs(
        scan.scan_dna(sample, [matrix], flat4, [threshold], window_size)[0]
    ) == baseline


def test_window_smaller_than_a_high_order_gram_is_rejected(flat4):
    matrix = [[0.0] * 3 for _ in range(16)]
    with pytest.raises(ValueError, match="q-gram"):
        scan.scan_dna("ACGTACGT", [matrix], flat4, [1.0], 1)


def test_window_size_must_be_positive(flat4, polya):
    with pytest.raises(ValueError, match="window size"):
        scan.Scanner(0)


def test_an_absurd_window_size_is_refused():
    with pytest.raises(ValueError, match="lookup table"):
        scan.Scanner(30)


# --------------------------------------------------------------------- Scanner
def test_scanner_matches_the_one_shot_function(sample, pfm_files, flat4):
    matrices = [parsers.pfm_to_log_odds(p, flat4, 0.01) for p in pfm_files]
    thresholds = [tools.threshold_from_p(m, flat4, 1e-3) for m in matrices]

    scanner = scan.Scanner(7)
    scanner.set_motifs(matrices, flat4, thresholds)
    assert len(scanner) == scanner.size() == len(matrices)

    for a, b in zip(
        scanner.scan(sample),
        scan.scan_dna(sample, matrices, flat4, thresholds),
        strict=True,
    ):
        assert as_pairs(a) == as_pairs(b)


def test_scanner_is_reusable_across_sequences(dna, pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-3)
    scanner = scan.Scanner(7)
    scanner.set_motifs([matrix], flat4, [threshold])

    for start in (0, 5000, 10000):
        chunk = dna[start : start + 5000]
        assert as_pairs(scanner.scan(chunk)[0]) == as_pairs(
            scan.scan_dna(chunk, [matrix], flat4, [threshold])[0]
        )


def test_an_empty_scanner_finds_nothing():
    assert scan.Scanner(7).scan("ACGTACGTACGT") == []
    assert len(scan.Scanner(7)) == 0


def test_set_motifs_can_be_called_again(sample, pfm_files, flat4):
    scanner = scan.Scanner(7)
    first = parsers.pfm_to_log_odds(pfm_files[0], flat4, 0.01)
    second = parsers.pfm_to_log_odds(pfm_files[1], flat4, 0.01)

    scanner.set_motifs([first], flat4, [tools.threshold_from_p(first, flat4, 1e-3)])
    scanner.set_motifs([second], flat4, [tools.threshold_from_p(second, flat4, 1e-3)])

    assert len(scanner) == 1
    expected = scan.scan_dna(
        sample, [second], flat4, [tools.threshold_from_p(second, flat4, 1e-3)]
    )
    assert as_pairs(scanner.scan(sample)[0]) == as_pairs(expected[0])


def test_background_does_not_change_which_positions_match(sample, pfm_file, flat4):
    """The scanner's bg only tunes the search; it must not affect the result."""
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-3)

    results = []
    for bg in (flat4, [0.4, 0.1, 0.1, 0.4], [0.1, 0.4, 0.4, 0.1]):
        scanner = scan.Scanner(7)
        scanner.set_motifs([matrix], bg, [threshold])
        results.append(as_pairs(scanner.scan(sample)[0]))
    assert results[0] == results[1] == results[2]


# -------------------------------------------------------------------- max hits
def test_scan_max_hits_keeps_the_leftmost(sample, pfm_file, flat4):
    """The cut-off stops the scan recording more, so it keeps the first ones
    found -- which is to say the leftmost, not the best-scoring."""
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-2)
    scanner = scan.Scanner(7)
    scanner.set_motifs([matrix], flat4, [threshold])

    everything = scanner.scan(sample)[0]
    limited = scanner.scan_max_hits(sample, 10)[0]
    assert len(everything) > 10
    assert as_pairs(limited) == as_pairs(everything[:10])


def test_counts_max_hits_agrees_with_scan_max_hits(sample, pfm_files, flat4):
    matrices = [parsers.pfm_to_log_odds(p, flat4, 0.01) for p in pfm_files[:4]]
    thresholds = [tools.threshold_from_p(m, flat4, 1e-2) for m in matrices]
    scanner = scan.Scanner(7)
    scanner.set_motifs(matrices, flat4, thresholds)
    assert scanner.counts_max_hits(sample, 25) == [
        len(h) for h in scanner.scan_max_hits(sample, 25)
    ]


def test_counts_without_a_limit_equals_scanning(sample, pfm_files, flat4):
    matrices = [parsers.pfm_to_log_odds(p, flat4, 0.01) for p in pfm_files[:4]]
    thresholds = [tools.threshold_from_p(m, flat4, 1e-2) for m in matrices]
    scanner = scan.Scanner(7)
    scanner.set_motifs(matrices, flat4, thresholds)
    huge = 10**9
    assert scanner.counts_max_hits(sample, huge) == [len(h) for h in scanner.scan(sample)]


def test_max_hits_must_not_be_negative(sample, pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    scanner = scan.Scanner(7)
    scanner.set_motifs([matrix], flat4, [0.0])
    with pytest.raises(ValueError, match="max_hits"):
        scanner.scan_max_hits("ACGT", -1)


# ------------------------------------------------------------------- best hits
def test_best_hits_returns_roughly_the_requested_number(dna, pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    target = 100
    hits = scan.scan_best_hits_dna(dna, [matrix], target)[0]
    assert target <= len(hits) <= 10 * target


def test_best_hits_are_the_highest_scoring_ones(dna, pfm_file, flat4):
    """Everything above the weakest returned hit must also have been returned."""
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    hits = scan.scan_best_hits_dna(dna, [matrix], 50)[0]
    cutoff = min(h.score for h in hits)
    exhaustive = scan.scan_dna(dna, [matrix], flat4, [cutoff])[0]
    assert as_pairs(hits) == as_pairs(exhaustive)


def test_best_hits_handles_several_matrices_at_once(dna, pfm_files, flat4):
    matrices = [parsers.pfm_to_log_odds(p, flat4, 0.01) for p in pfm_files]
    results = scan.scan_best_hits_dna(dna, matrices, 200)
    assert len(results) == len(matrices)


def test_best_hits_needs_a_positive_target(dna, pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    with pytest.raises(ValueError, match="target"):
        scan.scan_best_hits_dna(dna, [matrix], 0)


# --------------------------------------------------------------- naive_scan_dna
def test_naive_scan_matches_the_reference_implementation(sample, pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-3)
    assert as_pairs(scan.naive_scan_dna(sample, matrix, threshold)) == rounded(
        reference.brute_force_scan(sample, matrix, threshold)
    )


def test_naive_and_fast_scans_find_the_same_positions(sample, pfm_file, flat4):
    """Scores may differ in their last bits -- the two add up the columns in a
    different order -- but never which positions match."""
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4, 0.01)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-3)
    naive = scan.naive_scan_dna(sample, matrix, threshold)
    fast = scan.scan_dna(sample, [matrix], flat4, [threshold])[0]
    assert [m.pos for m in naive] == [m.pos for m in fast]
    for a, b in zip(naive, fast, strict=True):
        assert a.score == pytest.approx(b.score, abs=1e-12)


def test_naive_high_order_scan_matches_brute_force(sample, adm_file, flat4):
    matrix = parsers.adm_to_log_odds(adm_file, flat4, 0.01)
    threshold = tools.threshold_from_p(matrix, flat4, 1e-3, a=4)
    assert as_pairs(scan.naive_scan_dna(sample, matrix, threshold, 4)) == rounded(
        reference.brute_force_scan(sample, matrix, threshold)
    )


# --------------------------------------------------------------------- variants
def test_variant_matches_finds_hits_created_by_an_ambiguity_code(flat4, polya):
    seq = "GGAAAAWAAGG"
    #      01234567890
    # Only four consecutive As, and the W both breaks the run and splits the
    # sequence into two scannable regions -- so nothing matches as written.
    scanner = scan.Scanner(7)
    scanner.set_motifs([polya], flat4, [5.0])

    variants = tools.snp_variants(seq)
    assert scanner.scan(seq)[0] == []

    with_variants = scanner.variant_matches(seq, variants)[0]
    assert [h.pos for h in with_variants] == [2, 3, 4]
    for hit in with_variants:
        assert all(variants[i].modified_seq == "A" for i in hit.variants)


def test_variant_matches_excludes_hits_the_plain_sequence_already_has(flat4, polya):
    """A match that does not depend on a variant belongs to scan(), not here."""
    scanner = scan.Scanner(7)
    scanner.set_motifs([polya], flat4, [5.0])
    seq = "GGAAAAAWGG"
    assert [m.pos for m in scanner.scan(seq)[0]] == [2]
    assert all(h.variants for h in scanner.variant_matches(seq, tools.snp_variants(seq))[0])


def test_variant_matches_reports_indices_into_the_caller_s_list(flat4, polya):
    """Variants are sorted internally, so the indices reported have to be
    translated back to the order they were given in."""
    scanner = scan.Scanner(7)
    scanner.set_motifs([polya], flat4, [5.0])
    seq = "GGAAAAWAAGG"
    variants = list(reversed(tools.snp_variants(seq)))
    for hit in scanner.variant_matches(seq, variants)[0]:
        for i in hit.variants:
            assert variants[i].modified_seq == "A"


def test_variant_matches_is_empty_without_variants(flat4, polya):
    scanner = scan.Scanner(7)
    scanner.set_motifs([polya], flat4, [5.0])
    assert scanner.variant_matches("AAAAAAAAAA", []) == [[]]


def test_variant_matches_handles_a_deletion(flat4, polya):
    scanner = scan.Scanner(7)
    scanner.set_motifs([polya], flat4, [5.0])
    seq = "GGAAAGAAAGG"
    deletion = tools.Variant(5, 6, "")  # removes the G splitting the A runs
    assert scanner.variant_matches(seq, [deletion])[0], (
        "deleting the interrupting base should create a 6-long A run"
    )


def test_variant_matches_handles_an_insertion(flat4, polya):
    scanner = scan.Scanner(7)
    scanner.set_motifs([polya], flat4, [5.0])
    seq = "GGAAAAGG"
    insertion = tools.Variant(6, 6, "A")  # end == start inserts without replacing
    assert scanner.variant_matches(seq, [insertion])[0]


def test_variant_matches_ignores_deletions_at_the_sequence_edges(flat4, polya):
    scanner = scan.Scanner(7)
    scanner.set_motifs([polya], flat4, [5.0])
    seq = "GAAAAAG"
    assert scanner.variant_matches(seq, [tools.Variant(0, 1, "")]) == [[]]
    assert scanner.variant_matches(seq, [tools.Variant(6, 7, "")]) == [[]]


def test_max_depth_limits_how_many_variants_combine(flat4, polya):
    scanner = scan.Scanner(7)
    scanner.set_motifs([polya], flat4, [5.0])
    seq = "GGAAAWWGG"  # both Ws have to become A for the run to reach five
    variants = tools.snp_variants(seq)

    assert scanner.variant_matches(seq, variants, 1) == [[]]
    for depth in (2, 3):
        hits = scanner.variant_matches(seq, variants, depth)[0]
        assert [(h.pos, len(h.variants)) for h in hits] == [(2, 2)]


def test_max_depth_zero_means_no_limit(flat4, polya):
    scanner = scan.Scanner(7)
    scanner.set_motifs([polya], flat4, [5.0])
    seq = "GGAAAWWGG"
    variants = tools.snp_variants(seq)
    assert scanner.variant_matches(seq, variants, 0) == scanner.variant_matches(
        seq, variants, 2
    )


# ----------------------------------------------------------- custom alphabets
def test_custom_alphabet_scan():
    seq = "ACGTMMMMMACGT"
    matrix = [
        [0.0, 0.0, 0.0],  # A
        [0.0, 0.0, 0.0],  # C
        [0.0, 0.0, 0.0],  # G
        [0.0, 0.0, 0.0],  # T
        [1.0, 1.0, 1.0],  # M
    ]
    alphabet = ["aA", "cC", "gG", "tT", "mM"]
    hits = scan.scan(seq, [matrix], tools.flat_bg(5), [3.0], 3, alphabet)[0]
    assert [h.pos for h in hits] == [4, 5, 6]


def test_custom_alphabet_agrees_with_brute_force():
    """An alphabet of five is not a power of two, so some window codes name
    symbols that do not exist and have to be rejected."""
    seq = "MACGTMMACGTMMMACGTM"
    matrix = [[0.5, -1.0], [0.1, 0.2], [-0.3, 0.4], [0.7, -0.2], [1.0, 1.0]]
    threshold = 0.5
    hits = scan.scan(
        seq, [matrix], tools.flat_bg(5), [threshold], 4, ["aA", "cC", "gG", "tT", "mM"]
    )[0]
    expected = reference.brute_force_scan(seq, matrix, threshold, alphabet="ACGTM")
    assert as_pairs(hits) == rounded(expected)


# ------------------------------------------------------------------ validation
def test_threshold_count_must_match_matrix_count(flat4):
    matrix = [[1.0] * 4, [0.0] * 4, [0.0] * 4, [0.0] * 4]
    with pytest.raises(ValueError, match="2 matrices but 1 thresholds"):
        scan.scan_dna("ACGT", [matrix, matrix], flat4, [1.0])
    with pytest.raises(ValueError, match="thresholds"):
        scan.Scanner(7).set_motifs([matrix], flat4, [1.0, 2.0])


def test_ragged_matrix_is_rejected(flat4):
    with pytest.raises(ValueError, match="differing lengths"):
        scan.scan_dna("ACGT", [[[1.0, 2.0], [3.0]]], flat4, [1.0])


def test_non_ascii_sequence_is_rejected(flat4):
    matrix = [[1.0] * 4, [0.0] * 4, [0.0] * 4, [0.0] * 4]
    with pytest.raises(ValueError, match="ASCII"):
        scan.scan_dna("ACGTéACGT", [matrix], flat4, [1.0])


def test_match_is_iterable_and_comparable():
    match = scan.Match(3, 1.5)
    assert tuple(match) == (3, 1.5)
    assert len(match) == 2
    assert match == scan.Match(3, 1.5)


def test_variants_sort_by_position():
    assert tools.Variant(1, 2, "A") < tools.Variant(3, 4, "C")
