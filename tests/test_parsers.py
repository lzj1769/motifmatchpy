"""Matrix file parsing, including the failure modes."""

from __future__ import annotations

import gzip

import pytest

import reference
from motifmatchpy import parsers, tools
from motifmatchpy.parsers import ParseError


def test_pfm_has_one_row_per_base(pfm_file):
    matrix = parsers.pfm(pfm_file)
    assert len(matrix) == 4
    assert len({len(row) for row in matrix}) == 1


def test_pfm_reads_the_counts_verbatim(tmp_path):
    path = tmp_path / "toy.pfm"
    path.write_text("1 2 3\n4 5 6\n7 8 9\n10 11 12\n")
    assert parsers.pfm(path) == [
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
        [10.0, 11.0, 12.0],
    ]


def test_pfm_ignores_blank_lines(tmp_path):
    path = tmp_path / "gappy.pfm"
    path.write_text("1 2\n\n3 4\n\n\n5 6\n7 8\n\n")
    assert parsers.pfm(path) == [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]]


def test_pfm_to_log_odds_equals_transforming_the_counts(pfm_file, flat4):
    reference.assert_matrix_close(
        parsers.pfm_to_log_odds(pfm_file, flat4, 0.01),
        tools.log_odds(parsers.pfm(pfm_file), flat4, 0.01),
    )


def test_pfm_to_log_odds_honours_the_log_base(pfm_file, flat4):
    reference.assert_matrix_close(
        parsers.pfm_to_log_odds(pfm_file, flat4, 0.01, log_base=2),
        tools.log_odds(parsers.pfm(pfm_file), flat4, 0.01, log_base=2),
    )


def test_pfm_accepts_str_and_path(pfm_file, flat4):
    assert parsers.pfm(str(pfm_file)) == parsers.pfm(pfm_file)


# ---------------------------------------------------------------------- adm
def test_adm_terms_have_the_expected_shapes(adm_file):
    first_order = parsers.adm_1o_terms(adm_file)
    zero_order = parsers.adm_0o_terms(adm_file)
    assert len(first_order) == 16
    assert len(zero_order) == 4
    assert len({len(row) for row in first_order}) == 1


def test_adm_to_log_odds_equals_transforming_the_terms(adm_file, flat4):
    first_order = parsers.adm_1o_terms(adm_file)
    zero_order = parsers.adm_0o_terms(adm_file)
    low = [[row[0] for row in zero_order]]
    reference.assert_matrix_close(
        parsers.adm_to_log_odds(adm_file, flat4, 0.01),
        tools.log_odds_high_order(first_order, low, flat4, 0.01, 4),
    )


def test_adm_to_log_odds_produces_a_first_order_matrix(adm_file, flat4):
    matrix = parsers.adm_to_log_odds(adm_file, flat4, 0.01)
    assert len(matrix) == 16
    assert reference.q_gram_size(len(matrix), 4) == 2


# ------------------------------------------------------------------- errors
def test_missing_file_raises_file_not_found(tmp_path, flat4):
    with pytest.raises(FileNotFoundError, match="no such file"):
        parsers.pfm(tmp_path / "absent.pfm")
    with pytest.raises(FileNotFoundError):
        parsers.pfm_to_log_odds(tmp_path / "absent.pfm", flat4, 0.01)


def test_directory_raises_is_a_directory(tmp_path):
    with pytest.raises(IsADirectoryError):
        parsers.pfm(tmp_path)


def test_empty_file_is_a_parse_error(tmp_path):
    """The C++ implementation reads row 0 before checking the table is
    non-empty, so this input crashes it outright."""
    path = tmp_path / "empty.pfm"
    path.write_text("")
    with pytest.raises(ParseError, match="could not parse"):
        parsers.pfm(path)


def test_whitespace_only_file_is_a_parse_error(tmp_path):
    path = tmp_path / "blank.pfm"
    path.write_text("\n\n   \n\t\n")
    with pytest.raises(ParseError, match="could not parse"):
        parsers.pfm(path)


def test_non_numeric_file_is_a_parse_error(tmp_path):
    path = tmp_path / "prose.pfm"
    path.write_text("this is not a matrix\nnor is this\n")
    with pytest.raises(ParseError, match="could not parse"):
        parsers.pfm(path)


def test_ragged_matrix_is_a_parse_error(tmp_path):
    path = tmp_path / "ragged.pfm"
    path.write_text("1 2 3\n4 5\n6 7 8\n9 10 11\n")
    with pytest.raises(ParseError, match="could not parse"):
        parsers.pfm(path)


def test_gzipped_numeric_matrix_is_decompressed(tmp_path):
    """Compression is detected by content, just as for motif collections."""
    path = tmp_path / "compressed.pfm.gz"
    path.write_bytes(gzip.compress(b"1 2\n3 4\n5 6\n7 8\n"))
    assert parsers.pfm(path) == [[1, 2], [3, 4], [5, 6], [7, 8]]


def test_adm_with_the_wrong_row_count_is_a_parse_error(pfm_file):
    with pytest.raises(ParseError, match="could not parse"):
        parsers.adm_1o_terms(pfm_file)


@pytest.mark.parametrize("log_base", [0.5, 1.0])
def test_bad_log_base_is_rejected(pfm_file, flat4, log_base):
    with pytest.raises(ValueError, match="log_base"):
        parsers.pfm_to_log_odds(pfm_file, flat4, 0.01, log_base)


def test_bad_background_is_rejected(pfm_file):
    with pytest.raises(ValueError, match="4 entries"):
        parsers.pfm_to_log_odds(pfm_file, [0.5, 0.5], 0.01)


def test_parse_error_is_a_value_error(tmp_path):
    """So callers can catch it without importing the parsers module."""
    path = tmp_path / "empty.pfm"
    path.write_text("")
    assert issubclass(ParseError, ValueError)
    with pytest.raises(ValueError):
        parsers.pfm(path)


def test_read_table_stops_a_row_at_the_first_non_number(tmp_path):
    path = tmp_path / "mixed.txt"
    path.write_text("1 2 x 3\n4 5 6\n")
    assert parsers.read_table(path) == [[1.0, 2.0], [4.0, 5.0, 6.0]]


def test_read_table_drops_rows_with_no_numbers(tmp_path):
    path = tmp_path / "commented.txt"
    path.write_text("# a comment\n1 2\n\n3 4\n")
    assert parsers.read_table(path) == [[1.0, 2.0], [3.0, 4.0]]


def test_adm_rejects_a_nonsensical_alphabet_size(adm_file):
    with pytest.raises(ValueError, match="alphabet size"):
        parsers.adm_1o_terms(adm_file, 0)


def test_adm_zero_order_rows_may_be_shorter(tmp_path):
    """Only the first value of each zero-order row is used, so the rest of the
    row need not be there at all."""
    rows = ["1 2 3"] * 16 + ["7"] * 4
    path = tmp_path / "short_tail.adm"
    path.write_text("\n".join(rows) + "\n")
    assert len(parsers.adm_1o_terms(path)) == 16
    assert parsers.adm_0o_terms(path) == [[7.0], [7.0], [7.0], [7.0]]
