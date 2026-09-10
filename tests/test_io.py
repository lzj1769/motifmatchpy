"""Reading sequence files."""

from __future__ import annotations

import gzip

import pytest

from motifmatchpy import io


def write(tmp_path, name, text, compress=False):
    path = tmp_path / name
    if compress:
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)
    return path


def test_read_fasta_single_record(tmp_path):
    path = write(tmp_path, "one.fa", ">chr1\nACGT\nACGT\n")
    assert list(io.read_fasta(path)) == [("chr1", "ACGTACGT")]


def test_read_fasta_several_records(tmp_path):
    path = write(tmp_path, "many.fa", ">a\nAC\nGT\n>b\nTTTT\n>c\nG\n")
    assert list(io.read_fasta(path)) == [("a", "ACGT"), ("b", "TTTT"), ("c", "G")]


def test_read_fasta_name_stops_at_the_first_space(tmp_path):
    """The identifier is the first token; the rest of the header is description."""
    path = write(tmp_path, "desc.fa", ">chr1 5000:55000 human\nACGT\n")
    assert list(io.read_fasta(path)) == [("chr1", "ACGT")]


def test_read_fasta_can_keep_the_whole_header(tmp_path):
    path = write(tmp_path, "desc.fa", ">chr1 5000:55000\nACGT\n")
    assert list(io.read_fasta(path, full_header=True)) == [("chr1 5000:55000", "ACGT")]


def test_read_fasta_handles_an_empty_record(tmp_path):
    path = write(tmp_path, "empty_record.fa", ">a\n>b\nACGT\n")
    assert list(io.read_fasta(path)) == [("a", ""), ("b", "ACGT")]


def test_read_fasta_of_an_empty_file(tmp_path):
    assert list(io.read_fasta(write(tmp_path, "empty.fa", ""))) == []


def test_read_sequences_accepts_plain_text(tmp_path):
    path = write(tmp_path, "plain.txt", "ACGT\nACGT\n")
    assert list(io.read_sequences(path)) == [("plain.txt", "ACGTACGT")]


def test_read_sequences_detects_fasta(tmp_path):
    path = write(tmp_path, "seq.txt", ">x\nACGT\n")
    assert list(io.read_sequences(path)) == [("x", "ACGT")]


def test_read_sequences_skips_leading_blank_lines(tmp_path):
    path = write(tmp_path, "padded.fa", "\n\n>x\nACGT\n")
    assert list(io.read_sequences(path)) == [("x", "ACGT")]


def test_read_sequences_of_an_empty_file(tmp_path):
    assert list(io.read_sequences(write(tmp_path, "empty.fa", ""))) == []


def test_gzip_is_detected_from_the_content_not_the_name(tmp_path):
    """A compressed file named .fa still reads, and vice versa."""
    plain = write(tmp_path, "a.fa", ">x\nACGT\n")
    zipped = write(tmp_path, "b.fa", ">x\nACGT\n", compress=True)
    assert list(io.read_fasta(zipped)) == list(io.read_fasta(plain))


def test_gzip_with_the_usual_suffix(tmp_path):
    path = write(tmp_path, "c.fa.gz", ">x\nACGTACGT\n", compress=True)
    assert list(io.read_sequences(path)) == [("x", "ACGTACGT")]


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        list(io.read_fasta(tmp_path / "absent.fa"))


def test_the_example_file_reads(sequence_file):
    records = list(io.read_fasta(sequence_file))
    assert len(records) == 1
    name, seq = records[0]
    assert name == "chr1"
    assert len(seq) == 500000
