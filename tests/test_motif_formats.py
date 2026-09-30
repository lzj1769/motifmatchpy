"""Content-based collection detection and equivalent matrix conversion."""

import gzip

import pytest

from motifmatchpy import Motif, cli, parsers, read_motifs, tools

COUNTS = [[10, 0, 0], [0, 10, 0], [0, 0, 10], [0, 0, 0]]
JASPAR = ">MA0000.1 TEST\nA [10 0 0]\nC [0 10 0]\nG [0 0 10]\nT [0 0 0]\n"
MEME = """MEME version 5

ALPHABET= ACGT
strands: + -
Background letter frequencies
A 0.4 C 0.1 G 0.1 T 0.4

MOTIF MA0000.1 TEST
letter-probability matrix: alength= 4 w= 3 nsites= 10 E= 0
1 0 0 0
0 1 0 0
0 0 1 0
"""
TRANSFAC = """VV EXAMPLE
XX
//
AC MA0000.1
ID TEST
DE test motif
P0 A C G T
01 10 0 0 0 A
02 0 10 0 0 C
03 0 0 10 0 G
XX
//
"""


@pytest.mark.parametrize("data, kind", [(JASPAR, "jaspar"), (MEME, "meme"),
                                        (TRANSFAC, "transfac")])
@pytest.mark.parametrize("compressed", [False, True])
def test_auto_detection_ignores_suffix_and_preserves_scores(tmp_path, data, kind, compressed):
    path = tmp_path / "mislabelled.pfm"
    content = ("# comment\n\n" + data).encode()
    path.write_bytes(gzip.compress(content) if compressed else content)
    assert parsers.detect_format(path) == kind
    assert getattr(parsers, kind)(path) == [("MA0000.1_TEST", COUNTS)]
    for bg in ([0.25] * 4, [0.3, 0.2, 0.2, 0.3]):
        for base in (None, 2):
            motif = read_motifs([path], bg=bg, log_base=base)[0]
            expected = tuple(map(tuple, tools.log_odds(COUNTS, bg, 0.01, base)))
            assert motif.matrix == expected
            assert Motif.from_file(path).name == "MA0000.1_TEST"


@pytest.mark.parametrize("data", [JASPAR, MEME, TRANSFAC])
def test_collections_expand_in_api_and_all_cli_commands(tmp_path, capsys, data):
    path = tmp_path / "motifs.data"
    # Repeating the format header is unnecessary; retain only the second record.
    second = {JASPAR: JASPAR, MEME: MEME[MEME.index("MOTIF "):],
              TRANSFAC: TRANSFAC[TRANSFAC.index("AC "):]}[data]
    path.write_text(data + "\n" + second.replace("MA0000.1", "MA0001.1"))
    assert len(read_motifs([path])) == 2
    with pytest.raises(ValueError, match="use read_motifs"):
        Motif.from_file(path)
    with pytest.raises(ValueError, match="one motif per file"):
        read_motifs([path], names=["ambiguous"])
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nTTACGTT\n")
    regions = tmp_path / "regions.bed"
    regions.write_text("chr1 2 6\n")
    for command, args in [
        ("info", []), ("threshold", ["-p", "0.02"]),
        ("scan", ["--regions", str(regions), "--genome", str(genome), "-p", "0.02"]),
    ]:
        assert cli.main([command, "--matrices", str(path), *args]) == 0
        output = capsys.readouterr().out
        assert "MA0000.1_TEST" in output and "MA0001.1_TEST" in output
        if command == "scan":
            assert len(output.splitlines()) == 4
            assert all(row.split("\t")[:3] in [["chr1", "2", "5"], ["chr1", "3", "6"]]
                       for row in output.splitlines())


def test_meme_defaults_and_inferred_width(tmp_path):
    path = tmp_path / "minimal"
    path.write_text(MEME.replace("alength= 4 w= 3 nsites= 10 E= 0", ""))
    assert parsers.meme(path)[0][1] == [[v * 2 for v in row] for row in COUNTS]


def test_transfac_column_order_and_id_only(tmp_path):
    path = tmp_path / "tf"
    path.write_text("ID TF\nPO T G A C\n01 4 3 1 2 N\n//\n")
    assert parsers.transfac(path) == [("TF", [[1], [2], [3], [4]])]


@pytest.mark.parametrize("data", [
    MEME.replace("ACGT", "ACGU"), MEME.replace("alength= 4", "alength= 20"),
    MEME.replace("w= 3", "w= 4"), MEME.replace("nsites= 10", "nsites= 0"),
    MEME.replace("1 0 0 0", "0.1 0 0 0"), MEME.replace("1 0 0 0", "nan 0 0 0"),
    MEME.replace("1 0 0 0", "-1 2 0 0"),
    "MEME version 5\nMOTIF X\nlog-odds matrix: alength= 4 w= 1\n1 2 3 4\n",
    TRANSFAC.replace("P0 A C G T", "P0 A C G G"),
    TRANSFAC.replace("02 0 10", "04 0 10"), TRANSFAC.rstrip()[:-2],
    TRANSFAC.replace("01 10", "01 -10"), "unknown format\n1 2 3 4\n",
])
def test_malformed_files_fail_without_numeric_fallback(tmp_path, data):
    path = tmp_path / "bad.pfm"
    path.write_text(data)
    with pytest.raises(ValueError):
        read_motifs([path])


@pytest.mark.parametrize("data", [JASPAR, MEME, TRANSFAC])
def test_score_input_does_not_silently_read_counts(tmp_path, data):
    path = tmp_path / "scores"
    path.write_text(data)
    with pytest.raises(ValueError, match="use -m, not -S"):
        read_motifs([path], log_odds=False)


def test_numeric_formats_detected_by_shape_even_with_wrong_suffix(tmp_path, adm_file):
    pfm = tmp_path / "plain.meme"
    pfm.write_text("10 0 0\n0 10 0\n0 0 10\n0 0 0\n")
    adm = tmp_path / "model.pfm.gz"
    adm.write_bytes(gzip.compress(adm_file.read_bytes()))
    assert parsers.detect_format(pfm) == "numeric"
    assert read_motifs([pfm])[0].matrix == tuple(map(tuple, tools.log_odds(COUNTS, [0.25]*4, .01)))
    assert read_motifs([adm])[0].matrix == Motif.from_file(adm_file).matrix
