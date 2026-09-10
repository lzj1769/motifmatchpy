"""Matrix and background utilities, checked against reference.py."""

from __future__ import annotations

import math

import pytest

import reference
from motifmatchpy import parsers, tools


# ----------------------------------------------------------------- backgrounds
def test_flat_bg_is_uniform():
    assert tools.flat_bg(4) == [0.25] * 4
    assert tools.flat_bg(5) == pytest.approx([0.2] * 5)


@pytest.mark.parametrize("alphabet_size", [0, -1])
def test_flat_bg_rejects_empty_alphabet(alphabet_size):
    with pytest.raises(ValueError, match="alphabet_size"):
        tools.flat_bg(alphabet_size)


@pytest.mark.parametrize("ps", [0.0, 1.0, 10.0])
def test_bg_from_sequence_matches_reference(dna, ps):
    assert tools.bg_from_sequence_dna(dna, ps) == pytest.approx(
        reference.background(dna, ps)
    )


def test_bg_from_sequence_ignores_non_acgt():
    with_n = tools.bg_from_sequence_dna("ACGTNNNNNNACGT", 0.0)
    without_n = tools.bg_from_sequence_dna("ACGTACGT", 0.0)
    assert with_n == pytest.approx(without_n)
    assert sum(with_n) == pytest.approx(1.0)


def test_bg_from_sequence_accepts_bytes(dna):
    assert tools.bg_from_sequence_dna(dna.encode()) == pytest.approx(
        tools.bg_from_sequence_dna(dna)
    )


def test_bg_from_sequence_rejects_non_ascii():
    with pytest.raises(ValueError, match="ASCII"):
        tools.bg_from_sequence_dna("ACGTé")


# ---------------------------------------------------------------- transforms
def test_log_odds_matches_reference(pfm_file, flat4):
    counts = parsers.pfm(pfm_file)
    reference.assert_matrix_close(
        tools.log_odds(counts, flat4, 0.01), reference.log_odds(counts, flat4, 0.01)
    )


def test_log_odds_with_base_rescales(pfm_file, flat4):
    counts = parsers.pfm(pfm_file)
    natural = tools.log_odds(counts, flat4, 0.01)
    base2 = tools.log_odds(counts, flat4, 0.01, log_base=2)
    for row_n, row_2 in zip(natural, base2):
        assert row_2 == pytest.approx([v / math.log(2) for v in row_n])


def test_log_odds_uses_background(pfm_file):
    counts = parsers.pfm(pfm_file)
    skewed = [0.4, 0.1, 0.1, 0.4]
    reference.assert_matrix_close(
        tools.log_odds(counts, skewed, 0.01), reference.log_odds(counts, skewed, 0.01)
    )


@pytest.mark.parametrize("log_base", [0.5, 1.0])
def test_log_odds_rejects_bad_base(pfm_file, flat4, log_base):
    with pytest.raises(ValueError, match="log_base"):
        tools.log_odds(parsers.pfm(pfm_file), flat4, 0.01, log_base)


def test_log_odds_high_order_dispatches_like_moods(adm_file, flat4):
    """MOODS overloads log_odds; the 5-argument call must reach the high-order form."""
    first_order = parsers.adm_1o_terms(adm_file)
    zero_order = parsers.adm_0o_terms(adm_file)
    low = [[row[0] for row in zero_order]]
    via_overload = tools.log_odds(first_order, low, flat4, 0.01, 4)
    via_explicit = tools.log_odds_high_order(first_order, low, flat4, 0.01, 4)
    assert via_overload == via_explicit
    assert len(via_overload) == 16


def test_reverse_complement_matches_reference(pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4)
    reference.assert_matrix_close(
        tools.reverse_complement(matrix), reference.reverse_complement(matrix)
    )


def test_reverse_complement_with_alphabet_size_agrees_for_zero_order(pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4)
    reference.assert_matrix_close(
        tools.reverse_complement(matrix, 4), tools.reverse_complement(matrix)
    )


def test_reverse_complement_high_order_matches_reference(adm_file, flat4):
    matrix = parsers.adm_to_log_odds(adm_file, flat4)
    reference.assert_matrix_close(
        tools.reverse_complement(matrix, 4), reference.reverse_complement(matrix, 4)
    )


@pytest.mark.parametrize("name", ["MA0001.pfm", "MA0007.pfm"])
def test_reverse_complement_is_an_involution(matrix_dir, flat4, name):
    matrix = parsers.pfm_to_log_odds(matrix_dir / name, flat4)
    twice = tools.reverse_complement(tools.reverse_complement(matrix, 4), 4)
    reference.assert_matrix_close(twice, matrix)


def test_reverse_complement_high_order_is_an_involution(adm_file, flat4):
    matrix = parsers.adm_to_log_odds(adm_file, flat4)
    twice = tools.reverse_complement(tools.reverse_complement(matrix, 4), 4)
    reference.assert_matrix_close(twice, matrix)


# -------------------------------------------------------------- score bounds
def test_score_bounds_match_exhaustive_search(flat4):
    # 5 columns over 4 symbols is 1024 sequences: small enough to enumerate.
    counts = [
        [10, 1, 1, 1, 5],
        [1, 10, 1, 1, 5],
        [1, 1, 10, 1, 5],
        [1, 1, 1, 10, 5],
    ]
    matrix = tools.log_odds(counts, flat4, 0.01)
    assert tools.max_score(matrix) == pytest.approx(reference.max_score(matrix))
    assert tools.min_score(matrix) == pytest.approx(reference.min_score(matrix))


def test_score_bounds_with_alphabet_size_agree_for_zero_order(pfm_files, flat4):
    for path in pfm_files:
        matrix = parsers.pfm_to_log_odds(path, flat4)
        assert tools.max_score(matrix, 4) == pytest.approx(tools.max_score(matrix))
        assert tools.min_score(matrix, 4) == pytest.approx(tools.min_score(matrix))


def test_high_order_score_bounds_match_exhaustive_search(flat4, rng):
    # 16 rows x 3 columns is a 4-long motif: 256 sequences.
    matrix = [[rng.uniform(-3, 3) for _ in range(3)] for _ in range(16)]
    assert tools.max_score(matrix, 4) == pytest.approx(reference.max_score(matrix, 4))
    assert tools.min_score(matrix, 4) == pytest.approx(reference.min_score(matrix, 4))


def test_min_delta_is_the_smallest_top_two_gap():
    """min_delta is the narrowest gap between the best and second-best symbol
    of any column -- not the narrowest gap between any two symbols."""
    matrix = [[0.0, 5.0], [1.0, 5.5], [3.0, 9.0], [10.0, 20.0]]
    # Column 0: best 10, runner-up 3 -> 7. Column 1: best 20, runner-up 9 -> 11.
    assert tools.min_delta(matrix) == pytest.approx(7.0)


def test_min_delta_ignores_tied_maxima():
    """Two symbols sharing the best score count as one, so the gap is measured
    against the next distinct value rather than reported as zero."""
    assert tools.min_delta([[1.0], [1.0], [0.0], [0.0]]) == pytest.approx(1.0)


def test_score_bounds_bracket_every_observed_score(dna, pfm_file, flat4):
    from motifmatchpy import scan

    matrix = parsers.pfm_to_log_odds(pfm_file, flat4)
    lo, hi = tools.min_score(matrix), tools.max_score(matrix)
    hits = scan.scan_dna(dna, [matrix], flat4, [lo])[0]
    assert hits, "a threshold at min_score should match everywhere"
    assert all(lo - 1e-9 <= h.score <= hi + 1e-9 for h in hits)


# ------------------------------------------------------------------ p-values
def test_threshold_lies_within_the_score_range(pfm_files, flat4):
    for path in pfm_files:
        matrix = parsers.pfm_to_log_odds(path, flat4)
        threshold = tools.threshold_from_p(matrix, flat4, 1e-4)
        assert tools.min_score(matrix) <= threshold <= tools.max_score(matrix)


def test_threshold_decreases_as_p_grows(pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4)
    thresholds = [tools.threshold_from_p(matrix, flat4, p) for p in (1e-6, 1e-4, 1e-2, 0.5)]
    assert thresholds == sorted(thresholds, reverse=True)


def test_threshold_with_alphabet_size_agrees_for_zero_order(pfm_files, flat4):
    """moods-dna passes a=4 for every matrix, so both forms must agree."""
    for path in pfm_files:
        matrix = parsers.pfm_to_log_odds(path, flat4)
        plain = tools.threshold_from_p(matrix, flat4, 1e-4)
        with_a = tools.threshold_from_p(matrix, flat4, 1e-4, a=4)
        assert with_a == pytest.approx(plain, abs=1e-6)


def test_higher_precision_converges(pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4)
    coarse = tools.threshold_from_p(matrix, flat4, 1e-4, precision=100)
    fine = tools.threshold_from_p(matrix, flat4, 1e-4, precision=20000)
    default = tools.threshold_from_p(matrix, flat4, 1e-4)
    assert abs(fine - default) <= abs(coarse - default) + 1e-9


def test_threshold_precision_alias_matches_keyword(pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4)
    assert tools.threshold_from_p_with_precision(matrix, flat4, 1e-4, 5000) == (
        tools.threshold_from_p(matrix, flat4, 1e-4, precision=5000)
    )


@pytest.mark.parametrize("p", [0.0, -0.1, 1.5])
def test_threshold_rejects_p_outside_range(pfm_file, flat4, p):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4)
    with pytest.raises(ValueError, match="p-value"):
        tools.threshold_from_p(matrix, flat4, p)


def test_threshold_rejects_bad_precision(pfm_file, flat4):
    matrix = parsers.pfm_to_log_odds(pfm_file, flat4)
    with pytest.raises(ValueError, match="precision"):
        tools.threshold_from_p(matrix, flat4, 1e-4, precision=0)


# ------------------------------------------------------------------ variants
def test_snp_variants_expands_iupac_codes():
    variants = tools.snp_variants("AWC")
    assert [(v.start_pos, v.end_pos, v.modified_seq) for v in variants] == [
        (1, 2, "A"),
        (1, 2, "T"),
    ]


def test_snp_variants_handles_every_ambiguity_code():
    expected = {
        "W": "AT", "S": "CG", "M": "AC", "K": "GT", "R": "AG", "Y": "CT",
        "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    }
    for code, alleles in expected.items():
        for symbol in (code, code.lower()):
            variants = tools.snp_variants(symbol)
            assert [v.modified_seq for v in variants] == list(alleles)


def test_snp_variants_ignores_unambiguous_bases():
    assert tools.snp_variants("ACGTNacgtn") == []


def test_variants_order_by_position():
    a = tools.Variant(1, 2, "A")
    b = tools.Variant(3, 4, "C")
    assert a < b
    assert not (b < a)


# ---------------------------------------------------------------- validation
def test_matrix_with_ragged_rows_is_rejected(flat4):
    with pytest.raises(ValueError, match="differing lengths"):
        tools.max_score([[1.0, 2.0], [3.0]])


def test_empty_matrix_is_rejected():
    with pytest.raises(ValueError, match="no rows"):
        tools.max_score([])


@pytest.mark.parametrize(
    "bg, message",
    [
        ([0.25, 0.25, 0.25], "4 entries"),
        ([0.5, 0.5, -0.1, 0.1], "non-negative"),
        ([0.0, 0.0, 0.0, 0.0], "positive value"),
    ],
)
def test_bad_background_is_rejected(pfm_file, bg, message):
    counts = parsers.pfm(pfm_file)
    with pytest.raises(ValueError, match=message):
        tools.log_odds(counts, bg, 0.01)
