"""The high-level API: motifs with names, hits with strands."""

from __future__ import annotations

import pytest

import motifmatchpy as mm
from motifmatchpy import parsers, tools


# ------------------------------------------------------------------- Motif
def test_motif_from_pfm_file(pfm_file):
    motif = mm.Motif.from_file(pfm_file)
    assert motif.name == "MA0001.pfm"
    assert motif.order == 0
    assert motif.length == motif.width == len(motif.matrix[0])
    assert len(motif.matrix) == 4
    assert len(motif) == motif.length


def test_motif_from_adm_file(adm_file):
    motif = mm.Motif.from_file(adm_file)
    assert motif.name == "E2F3.adm"
    assert motif.order == 1
    assert len(motif.matrix) == 16
    # A first-order model spans one position more than it has columns.
    assert motif.length == motif.width + 1


def test_motif_detects_the_format_regardless_of_suffix(tmp_path, adm_file):
    disguised = tmp_path / "mislabelled.pfm"
    disguised.write_text(adm_file.read_text())
    assert mm.Motif.from_file(disguised).order == 1


def test_motif_score_matrix_is_used_as_is(pfm_file):
    raw = mm.Motif.from_file(pfm_file, log_odds=False)
    assert [list(row) for row in raw.matrix] == parsers.pfm(pfm_file)


def test_motif_reverse_complement_flips_the_strand(pfm_file):
    motif = mm.Motif.from_file(pfm_file)
    rc = motif.reverse_complement()
    assert rc.strand == "-"
    assert rc.name == motif.name
    assert rc.reverse_complement().matrix == motif.matrix
    assert rc.max_score == pytest.approx(motif.max_score)


def test_motif_scores_match_the_tools(pfm_file, flat4):
    motif = mm.Motif.from_file(pfm_file)
    assert motif.max_score == tools.max_score(motif.matrix, 4)
    assert motif.min_score == tools.min_score(motif.matrix, 4)
    assert motif.threshold_from_p(1e-4) == tools.threshold_from_p(
        motif.matrix, flat4, 1e-4, a=4
    )


def test_motif_is_hashable_and_comparable(pfm_file):
    a = mm.Motif.from_file(pfm_file)
    b = mm.Motif.from_file(pfm_file)
    assert a == b
    assert len({a, b}) == 1


@pytest.mark.parametrize(
    "matrix, message",
    [
        ([], "no rows"),
        ([[], [], [], []], "no columns"),
        ([[1.0, 2.0], [3.0]], "differing lengths"),
        ([[1.0], [2.0], [3.0]], r"4\*\*k rows"),
    ],
)
def test_motif_rejects_bad_matrices(matrix, message):
    with pytest.raises(ValueError, match=message):
        mm.Motif(name="bad", matrix=matrix)


def test_motif_rejects_a_bad_strand():
    with pytest.raises(ValueError, match="strand"):
        mm.Motif(name="x", matrix=[[1.0], [1.0], [1.0], [1.0]], strand="?")


def test_read_motifs(pfm_files):
    motifs = mm.read_motifs(pfm_files)
    assert len(motifs) == len(pfm_files)
    assert [m.name for m in motifs] == [p.name for p in pfm_files]


def test_read_motifs_accepts_explicit_names(pfm_files):
    names = [f"motif{i}" for i in range(len(pfm_files))]
    assert [m.name for m in mm.read_motifs(pfm_files, names=names)] == names


def test_read_motifs_checks_the_name_count(pfm_files):
    with pytest.raises(ValueError, match="names"):
        mm.read_motifs(pfm_files, names=["only-one"])


def test_read_motifs_reports_an_unparseable_file(tmp_path):
    bad = tmp_path / "bad.pfm"
    bad.write_text("not a matrix at all\n")
    with pytest.raises(ValueError, match="could not parse"):
        mm.read_motifs([bad])


def test_read_motifs_reports_a_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        mm.read_motifs([tmp_path / "absent.pfm"])


# ----------------------------------------------------------- MotifScanner
def test_scanner_finds_hits_on_both_strands(sample, pfm_file):
    motifs = [mm.Motif.from_file(pfm_file)]
    both = mm.MotifScanner(motifs, p_value=1e-3)
    forward = mm.MotifScanner(motifs, p_value=1e-3, both_strands=False)

    both_hits = both.scan(sample)
    forward_hits = forward.scan(sample)
    assert {h.strand for h in both_hits} == {"+", "-"}
    assert {h.strand for h in forward_hits} == {"+"}
    assert len(both_hits) > len(forward_hits)


def test_scanner_hits_are_ordered_by_position(sample, pfm_files):
    scanner = mm.MotifScanner(mm.read_motifs(pfm_files), p_value=1e-4)
    hits = scanner.scan(sample)
    assert [h.pos for h in hits] == sorted(h.pos for h in hits)


def test_scanner_matches_the_low_level_scan(sample, pfm_file, flat4):
    from motifmatchpy import scan as low

    motif = mm.Motif.from_file(pfm_file)
    threshold = motif.threshold_from_p(1e-3)
    scanner = mm.MotifScanner([motif], p_value=1e-3, both_strands=False)
    hits = scanner.scan(sample)
    expected = low.scan_dna(sample, [motif.matrix], flat4, [threshold])[0]
    assert [(h.pos, h.score) for h in hits] == [(m.pos, m.score) for m in expected]


def test_scanner_with_an_absolute_threshold(sample, pfm_file):
    motif = mm.Motif.from_file(pfm_file)
    scanner = mm.MotifScanner([motif], threshold=8.0, both_strands=False)
    hits = scanner.scan(sample)
    assert hits
    assert all(h.score >= 8.0 for h in hits)
    assert scanner.thresholds() == [8.0]


def test_scanner_with_an_estimated_background(sample, pfm_file):
    """``bg='auto'`` re-derives the threshold from each sequence's own
    composition, so it must differ from the flat-background one."""
    motif = mm.Motif.from_file(pfm_file)
    auto = mm.MotifScanner([motif], p_value=1e-3, bg="auto")
    flat = mm.MotifScanner([motif], p_value=1e-3, bg="flat")
    assert auto.thresholds(sample) != flat.thresholds()
    assert auto.scan(sample)


def test_scanner_with_an_explicit_background(sample, pfm_file, skewed4):
    motif = mm.Motif.from_file(pfm_file)
    scanner = mm.MotifScanner([motif], p_value=1e-3, bg=skewed4)
    assert scanner.thresholds() == [motif.threshold_from_p(1e-3, skewed4)]


def test_scanner_counts_both_strands(sample, pfm_files):
    motifs = mm.read_motifs(pfm_files)
    scanner = mm.MotifScanner(motifs, p_value=1e-3)
    counts = scanner.count(sample)
    hits = scanner.scan(sample)
    assert len(counts) == len(motifs)
    assert sum(counts) == len(hits)


def test_scanner_max_hits_caps_each_motif(sample, pfm_files):
    scanner = mm.MotifScanner(mm.read_motifs(pfm_files), p_value=1e-2)
    limited = scanner.scan(sample, max_hits=5)
    per_motif: dict[tuple[str, str], int] = {}
    for hit in limited:
        per_motif[(hit.name, hit.strand)] = per_motif.get((hit.name, hit.strand), 0) + 1
    assert max(per_motif.values()) <= 5


def test_scanner_reports_variant_hits(flat4):
    motif = mm.Motif(name="polyA", matrix=[[1.0] * 5, [0.0] * 5, [0.0] * 5, [0.0] * 5])
    scanner = mm.MotifScanner([motif], threshold=5.0, both_strands=False)
    seq = "GGAAAAWAAGG"
    plain = scanner.scan(seq)
    with_variants = scanner.scan(seq, include_variants=True)
    assert not plain
    assert with_variants
    assert all(h.variants for h in with_variants)


def test_scan_best_hits_needs_no_threshold(dna, pfm_file):
    motif = mm.Motif.from_file(pfm_file)
    scanner = mm.MotifScanner([motif], both_strands=False)
    hits = scanner.scan_best_hits(dna, 100)
    assert 100 <= len(hits) <= 1000


def test_scan_without_a_threshold_is_refused(sample, pfm_file):
    scanner = mm.MotifScanner([mm.Motif.from_file(pfm_file)])
    with pytest.raises(ValueError, match="p_value or threshold"):
        scanner.scan(sample)


def test_scanner_rejects_two_cutoffs(pfm_file):
    with pytest.raises(ValueError, match="not both"):
        mm.MotifScanner([mm.Motif.from_file(pfm_file)], p_value=1e-3, threshold=5.0)


def test_scanner_rejects_an_empty_motif_list():
    with pytest.raises(ValueError, match="at least one motif"):
        mm.MotifScanner([])


def test_scanner_rejects_mixed_alphabets():
    a = mm.Motif(name="dna", matrix=[[1.0], [0.0], [0.0], [0.0]])
    b = mm.Motif(name="five", matrix=[[1.0], [0.0], [0.0], [0.0], [0.0]], alphabet_size=5)
    with pytest.raises(ValueError, match="alphabet size"):
        mm.MotifScanner([a, b], threshold=1.0)


def test_auto_background_needs_a_p_value(pfm_file):
    with pytest.raises(ValueError, match="bg='auto'"):
        mm.MotifScanner([mm.Motif.from_file(pfm_file)], threshold=5.0, bg="auto")


def test_auto_background_thresholds_need_a_sequence(pfm_file):
    scanner = mm.MotifScanner([mm.Motif.from_file(pfm_file)], p_value=1e-3, bg="auto")
    with pytest.raises(ValueError, match="needs a sequence"):
        scanner.thresholds()


def test_unknown_background_keyword_is_refused(pfm_file):
    with pytest.raises(ValueError, match="'flat' or 'auto'"):
        mm.MotifScanner([mm.Motif.from_file(pfm_file)], p_value=1e-3, bg="uniform")


# --------------------------------------------------------------------- Hit
def test_hit_reports_its_span_and_sequence(sample, pfm_file):
    motif = mm.Motif.from_file(pfm_file)
    scanner = mm.MotifScanner([motif], p_value=1e-3, both_strands=False)
    hit = scanner.scan(sample, sequence_name="chr1")[0]
    assert hit.name == motif.name
    assert hit.end == hit.pos + motif.length
    assert hit.matched_sequence(sample) == sample[hit.pos : hit.end]
    assert len(hit.matched_sequence(sample)) == motif.length
    assert hit.sequence_name == "chr1"


def test_hit_variant_sequence_spells_out_substitutions(flat4):
    motif = mm.Motif(name="polyA", matrix=[[1.0] * 5, [0.0] * 5, [0.0] * 5, [0.0] * 5])
    scanner = mm.MotifScanner([motif], threshold=5.0, both_strands=False)
    seq = "GGAAAAWAAGG"
    hit = scanner.scan(seq, include_variants=True)[0]
    rendered = hit.variant_sequence(seq)
    assert rendered.upper() == "AAAAA"
    # The substituted base is the only upper-case one.
    assert sum(c.isupper() for c in rendered) == len(hit.variants)


def test_hit_without_variants_renders_the_plain_window(sample, pfm_file):
    scanner = mm.MotifScanner([mm.Motif.from_file(pfm_file)], p_value=1e-3)
    hit = scanner.scan(sample)[0]
    assert hit.variants == ()
    assert hit.variant_sequence(sample) == hit.matched_sequence(sample).lower()


# ------------------------------------------------------------------ FASTA
def test_scan_fasta_streams_every_record(sequence_file, pfm_file):
    scanner = mm.MotifScanner([mm.Motif.from_file(pfm_file)], p_value=1e-4)
    hits = list(scanner.scan_fasta(sequence_file))
    assert hits
    assert {h.sequence_name for h in hits} == {"chr1"}


def test_scan_fasta_matches_scanning_the_string(sequence_file, dna, pfm_file):
    scanner = mm.MotifScanner([mm.Motif.from_file(pfm_file)], p_value=1e-4)
    from_file = list(scanner.scan_fasta(sequence_file))
    direct = scanner.scan(dna, sequence_name="chr1")
    assert [(h.name, h.pos, h.strand, h.score) for h in from_file] == [
        (h.name, h.pos, h.strand, h.score) for h in direct
    ]
