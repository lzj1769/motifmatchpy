"""Peak scanning: genomic coordinates and multi-motif JASPAR input."""

from __future__ import annotations

import gzip

import pytest

from motifmatchpy import cli, parsers
from motifmatchpy._regions import peak_sequences, read_peaks

JASPAR = """>MA0000.1 TEST
A [10 0 0]
C [0 10 0]
G [0 0 10]
T [0 0 0]
>MA0001.1 OTHER
0 0 10
0 10 0
10 0 0
0 0 0
"""


@pytest.fixture
def inputs(tmp_path):
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1 description\nTTACGT\nTTACGAA\n>chr2\nGGACGCC\n")
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t8\t11\tpeak2\t0\t-\nchr2\t1\t6\nchr1\t2\t6\n")
    scores = tmp_path / "scores.pfm"
    scores.write_text("2 -1 -1\n-1 2 -1\n-1 -1 2\n-1 -1 -1\n")
    jaspar = tmp_path / "motifs.jaspar"
    jaspar.write_text(JASPAR)
    return genome, peaks, scores, jaspar


def run(capsys, *args):
    status = cli.main([str(arg) for arg in args])
    capture = capsys.readouterr()
    return status, capture.out, capture.err


def test_genomic_coordinates_both_strands_and_peak_boundaries(inputs, capsys):
    genome, peaks, scores, _ = inputs
    status, out, err = run(
        capsys, "scan", "--regions", peaks, "--genome", genome, "-S", scores, "-t", 6,
    )
    assert status == 0, err
    assert out.splitlines() == [
        "chr1\t2\t5\tscores.pfm\t6\t+",
        "chr1\t3\t6\tscores.pfm\t6\t-",
        "chr1\t8\t11\tscores.pfm\t6\t+",
        "chr2\t2\t5\tscores.pfm\t6\t+",
    ]


def test_jaspar_collection_scan_and_other_subcommands(inputs, capsys):
    genome, peaks, _, motifs = inputs
    status, out, err = run(
        capsys, "scan", "--regions", peaks, "--genome", genome,
        "-m", motifs, "-p", 0.02,
    )
    assert status == 0, err
    assert len(out.splitlines()) == 4
    assert all(row.split("\t")[3] == "MA0000.1_TEST" for row in out.splitlines())
    for command, extra in [("info", []), ("threshold", ["-p", "0.02"])]:
        status, out, err = run(capsys, command, "-m", motifs, *extra)
        assert status == 0, err
        assert [row.split(",")[0] for row in out.splitlines()] == [
            "MA0000.1_TEST", "MA0001.1_OTHER",
        ]


def test_peak_default_pvalue_and_fixed_background(inputs, capsys):
    genome, peaks, scores, _ = inputs
    args = ["scan", "--regions", peaks, "--genome", genome, "-S", scores]
    default = run(capsys, *args)
    explicit = run(capsys, *args, "-p", "1e-4", "--batch")
    assert default == explicit
    assert default[0] == 0


def test_csv_offsets_leave_match_sequence_local(inputs, capsys):
    genome, peaks, scores, _ = inputs
    status, out, err = run(
        capsys, "scan", "--regions", peaks, "--genome", genome,
        "-S", scores, "-t", 6, "-R", "-f", "csv",
    )
    assert status == 0, err
    assert out.splitlines()[0] == "chr1,scores.pfm,2,+,6,ACG,"
    assert len(out.splitlines()) == 3


def test_overlapping_peaks_are_independent_and_never_joined(inputs, capsys):
    genome, peaks, scores, _ = inputs
    peaks.write_text("chr1 2 5\nchr1 2 5\nchr1 8 10\nchr1 10 11\n")
    status, out, err = run(
        capsys, "scan", "--regions", peaks, "--genome", genome, "-S", scores, "-t", 6,
    )
    assert status == 0, err
    assert out.splitlines() == ["chr1\t2\t5\tscores.pfm\t6\t+"] * 2


@pytest.mark.parametrize("peak, message", [
    ("chr1 -1 4", "expected BED"), ("chr1 4 4", "expected BED"),
    ("chr1 4", "expected BED"), ("chr1 x 4", "expected BED"),
    ("chr1 0 100", "exceeds reference length"), ("missing 0 3", "missing from reference"),
])
def test_bad_peak_errors(inputs, capsys, peak, message):
    genome, peaks, scores, _ = inputs
    peaks.write_text(peak + "\n")
    status, _, err = run(
        capsys, "scan", "--regions", peaks, "--genome", genome, "-S", scores, "-t", 6,
    )
    assert status == 2
    assert message in err


@pytest.mark.parametrize("options", [
    ["--regions", "peaks.bed"], ["--genome", "genome.fa"],
    ["--regions", "peaks.bed", "--genome", "genome.fa", "-s", "seq.fa"], [],
])
def test_conflicting_or_missing_sequence_inputs(capsys, options):
    status, _, _ = run(capsys, "scan", *options, "-t", 6)
    assert status == 2


def test_output_cannot_replace_peak_input(inputs, capsys):
    genome, peaks, scores, _ = inputs
    original = peaks.read_text()
    status, _, err = run(
        capsys, "scan", "--regions", peaks, "--genome", genome,
        "-S", scores, "-t", 6, "-o", peaks,
    )
    assert status == 2
    assert "overwrite" in err
    assert peaks.read_text() == original


def test_gzip_inputs_and_bed_metadata(inputs, tmp_path):
    genome, _, _, motifs = inputs
    peaks = tmp_path / "peaks.gz"
    with gzip.open(peaks, "wt") as handle:
        handle.write(
            "track name=test\nbrowser position chr1\n# comment\nchr1 2 5 peak 0 - 1 2 3 4\n"
        )
    compressed = tmp_path / "genome.gz"
    with gzip.open(compressed, "wt") as handle:
        handle.write(genome.read_text())
    records = list(peak_sequences(compressed, read_peaks(peaks)))
    assert [(region.start, seq) for region, seq in records] == [(2, "ACG")]
    with gzip.open(tmp_path / "motifs.gz", "wt") as handle:
        handle.write(motifs.read_text())
    assert parsers.jaspar(tmp_path / "motifs.gz") == parsers.jaspar(motifs)


@pytest.mark.parametrize("contents", [
    "", ">\n", "1 2 3\n", ">id\n1 2\n1\n1\n1\n",
    ">id\nA [1]\nG [1]\nC [1]\nT [1]\n",
    ">id\n1\n1\n1\n-1\n", ">id\n1\n1\n1\nnan\n",
    ">id\n1\n1\n1\n1\n>id\n1\n1\n1\n1\n",
])
def test_malformed_jaspar_rejected(tmp_path, contents):
    filename = tmp_path / "bad.jaspar"
    filename.write_text(contents)
    with pytest.raises(parsers.ParseError):
        parsers.jaspar(filename)


def test_empty_peak_list_has_no_hits(inputs, capsys):
    genome, peaks, scores, _ = inputs
    peaks.write_text("")
    status, out, err = run(
        capsys, "scan", "--regions", peaks, "--genome", genome, "-S", scores, "-t", 6,
    )
    assert status == 0, err
    assert out == ""


@pytest.mark.parametrize("reference, message", [
    (">chr1\nACG\n>chr1\nACG\n", "duplicate FASTA"),
    ("ACGT", "expected a FASTA header"),
    (">\nACG", "empty FASTA header"),
    (">chr1\nAC GT", "whitespace/ASCII"),
])
def test_invalid_reference(inputs, capsys, reference, message):
    genome, peaks, scores, _ = inputs
    genome.write_text(reference)
    peaks.write_text("chr1 0 3\n")
    status, _, err = run(
        capsys, "scan", "--regions", peaks, "--genome", genome, "-S", scores, "-t", 6,
    )
    assert status == 2
    assert message in err


def test_jaspar_also_works_in_original_sequence_mode(inputs, capsys):
    genome, _, _, motifs = inputs
    genome.write_text(">toy\nACGNNGCA\n")
    status, out, err = run(
        capsys, "scan", "-s", genome, "-m", motifs, "-p", 0.02, "--batch", "-R",
    )
    assert status == 0, err
    assert [line.split(",")[:4] for line in out.splitlines()] == [
        ["toy", "MA0000.1_TEST", "0", "+"], ["toy", "MA0001.1_OTHER", "5", "+"],
    ]
