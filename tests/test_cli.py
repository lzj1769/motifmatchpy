"""The command line interface."""

from __future__ import annotations

import csv
import io as _io

import pytest

from motifmatchpy import cli


def run(capsys, *args):
    """Run the CLI and return ``(exit_code, stdout, stderr)``."""
    code = cli.main([str(a) for a in args])
    captured = capsys.readouterr()
    return code, captured.out, captured.err


def rows(text, delimiter=","):
    return list(csv.reader(_io.StringIO(text), delimiter=delimiter))


# ----------------------------------------------------------------- plumbing
def test_no_command_prints_help(capsys):
    code, _, err = run(capsys)
    assert code == 1
    assert "usage" in err.lower()


def test_version(capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["--version"])
    assert exit_info.value.code == 0
    assert "motifmatchpy" in capsys.readouterr().out


# --------------------------------------------------------------------- info
def test_info(capsys, pfm_file, adm_file):
    code, out, _ = run(capsys, "info", "-m", pfm_file, adm_file, "--header")
    assert code == 0
    table = rows(out)
    assert table[0] == list(cli.INFO_COLUMNS)
    assert [r[0] for r in table[1:]] == ["MA0001.pfm", "E2F3.adm"]
    assert [r[2] for r in table[1:]] == ["0", "1"]  # model order


def test_info_without_a_header(capsys, pfm_file):
    _, out, _ = run(capsys, "info", "-m", pfm_file)
    assert rows(out)[0][0] == "MA0001.pfm"


# ---------------------------------------------------------------- threshold
def test_threshold(capsys, pfm_file):
    code, out, _ = run(capsys, "threshold", "-m", pfm_file, "-p", "1e-4", "--header")
    assert code == 0
    table = rows(out)
    assert table[0] == ["motif", "threshold", "max_score", "min_score"]
    assert float(table[1][1]) == pytest.approx(7.2695)


def test_threshold_rejects_a_bad_p_value(capsys, pfm_file):
    code, _, err = run(capsys, "threshold", "-m", pfm_file, "-p", "0")
    assert code == 2
    assert "-p must be in (0, 1]" in err


# --------------------------------------------------------------------- scan
def test_scan_csv(capsys, pfm_file, sequence_file):
    code, out, _ = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "1e-5", "--header"
    )
    assert code == 0
    table = rows(out)
    assert table[0] == list(cli.SCAN_COLUMNS)
    assert table[1:], "expected at least one match"
    for row in table[1:]:
        assert len(row) == len(cli.SCAN_COLUMNS)
        assert row[0] == "chr1"
        assert row[3] in ("+", "-")


def test_scan_tsv(capsys, pfm_file, sequence_file):
    _, out, _ = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "1e-5", "-f", "tsv"
    )
    assert "\t" in out.splitlines()[0]


def test_scan_bed(capsys, pfm_file, sequence_file):
    _, out, _ = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "1e-5", "-f", "bed"
    )
    for row in rows(out, "\t"):
        chrom, start, end, _name, score, strand = row
        assert chrom == "chr1"
        assert int(end) > int(start)
        assert strand in ("+", "-")
        float(score)


def test_bed_spans_match_the_motif_length(capsys, pfm_file, sequence_file):
    import motifmatchpy as mm

    length = mm.Motif.from_file(pfm_file).length
    _, out, _ = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "1e-5", "-f", "bed"
    )
    assert all(int(r[2]) - int(r[1]) == length for r in rows(out, "\t"))


def test_scan_with_an_absolute_threshold(capsys, pfm_file, sequence_file):
    _, out, _ = run(capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-t", "9")
    assert all(float(r[4]) >= 9 for r in rows(out))


def test_scan_forward_strand_only(capsys, pfm_file, sequence_file):
    _, out, _ = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "1e-5", "-R"
    )
    assert {r[3] for r in rows(out)} == {"+"}


def test_scan_best_hits(capsys, pfm_file, sequence_file):
    _, out, _ = run(capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-B", "20")
    assert 20 <= len(rows(out)) <= 400


def test_scan_batch_uses_the_given_background(capsys, pfm_file, sequence_file):
    """--batch keeps the thresholds fixed instead of re-deriving them per
    sequence, so a skewed --bg has to change the output."""
    _, flat, _ = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "1e-4", "--batch"
    )
    _, skewed, _ = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "1e-4", "--batch",
        "--bg", "0.4", "0.1", "0.1", "0.4",
    )
    assert flat != skewed


def test_scan_score_matrices_are_used_directly(capsys, pfm_file, sequence_file):
    """-S skips the log-odds conversion, so the same file gives different scores."""
    _, converted, _ = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-t", "9"
    )
    _, direct, _ = run(
        capsys, "scan", "-S", pfm_file, "-s", sequence_file, "-t", "9"
    )
    assert converted != direct


def test_scan_writes_to_a_file(capsys, tmp_path, pfm_file, sequence_file):
    out_path = tmp_path / "hits.csv"
    code, out, _ = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "1e-5",
        "-o", out_path,
    )
    assert code == 0
    assert out == ""
    assert len(rows(out_path.read_text())) > 0


def test_scan_reports_progress_when_asked(capsys, pfm_file, sequence_file):
    _, _, err = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "1e-5", "-vv"
    )
    assert "scanning chr1" in err
    assert "matches in chr1" in err


def test_scan_handles_several_sequence_files(capsys, tmp_path, pfm_file):
    first = tmp_path / "a.fa"
    first.write_text(">a\n" + "ACGTTTACGTAAATTTGGGCCC" * 20 + "\n")
    second = tmp_path / "b.fa"
    second.write_text(">b\n" + "TTTTAAAACCCCGGGG" * 20 + "\n")
    _, out, _ = run(
        capsys, "scan", "-m", pfm_file, "-s", first, second, "--threshold=-5"
    )
    assert {r[0] for r in rows(out)} == {"a", "b"}


def test_scan_reports_variant_matches(capsys, tmp_path, pfm_file):
    seq_path = tmp_path / "iupac.fa"
    # One ambiguity code per repeat, far enough apart that a 10 bp motif still
    # has room to span it.
    seq_path.write_text(">amb\n" + "ACGTAAATTTGGGCCCACGTWTTTGGGCCCAAATTT" * 10 + "\n")

    _, out, _ = run(capsys, "scan", "-m", pfm_file, "-s", seq_path, "-t", "2")
    assert any(r[6] for r in rows(out)), "expected some variant-dependent matches"

    _, without, _ = run(
        capsys, "scan", "-m", pfm_file, "-s", seq_path, "-t", "2", "--no-snps"
    )
    assert all(not r[6] for r in rows(without))
    assert len(rows(without)) < len(rows(out))


def test_fields_containing_the_separator_are_quoted(capsys, tmp_path, pfm_file):
    seq_path = tmp_path / "comma.fa"
    seq_path.write_text(">has,comma\n" + "ACGTAAATTTGGGCCC" * 20 + "\n")
    _, out, _ = run(capsys, "scan", "-m", pfm_file, "-s", seq_path, "-t", "0")
    assert all(r[0] == "has,comma" for r in rows(out))


# ------------------------------------------------------------------- errors
def test_scan_needs_a_cutoff(capsys, pfm_file, sequence_file):
    code, _, err = run(capsys, "scan", "-m", pfm_file, "-s", sequence_file)
    assert code == 2
    assert "no score cutoff" in err


def test_scan_refuses_two_cutoffs(capsys, pfm_file, sequence_file):
    code, _, err = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "0.001", "-t", "5"
    )
    assert code == 2
    assert "only one of" in err


def test_scan_needs_a_matrix(capsys, sequence_file):
    code, _, err = run(capsys, "scan", "-s", sequence_file, "-p", "0.001")
    assert code == 2
    assert "no matrix files" in err


def test_scan_reports_a_missing_matrix(capsys, tmp_path, sequence_file):
    code, _, err = run(
        capsys, "scan", "-m", tmp_path / "absent.pfm", "-s", sequence_file, "-p", "0.001"
    )
    assert code == 2
    assert "no such file" in err


def test_scan_reports_a_missing_sequence_file(capsys, tmp_path, pfm_file):
    code, _, err = run(
        capsys, "scan", "-m", pfm_file, "-s", tmp_path / "absent.fa", "-p", "0.001"
    )
    assert code == 2
    assert "could not read" in err


def test_scan_rejects_a_bad_log_base(capsys, pfm_file, sequence_file):
    code, _, err = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "0.001",
        "--log-base", "1",
    )
    assert code == 2
    assert "--log-base" in err


def test_scan_rejects_a_bad_p_value(capsys, pfm_file, sequence_file):
    code, _, err = run(capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "2")
    assert code == 2
    assert "-p must be in (0, 1]" in err


def test_scan_rejects_a_bad_best_hits_count(capsys, pfm_file, sequence_file):
    code, _, err = run(capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-B", "0")
    assert code == 2
    assert "-B must be at least 1" in err


def test_unwritable_output_is_reported(capsys, tmp_path, pfm_file, sequence_file):
    code, _, err = run(
        capsys, "scan", "-m", pfm_file, "-s", sequence_file, "-p", "0.001",
        "-o", tmp_path,  # a directory
    )
    assert code == 2
    assert "could not open" in err
