"""The bit-packing primitives."""

from __future__ import annotations

import pytest

from motifmatchpy import misc
from motifmatchpy._alphabet import build_map


@pytest.mark.parametrize(
    "alphabet_size, expected", [(1, 0), (2, 1), (3, 2), (4, 2), (5, 3), (8, 3), (9, 4)]
)
def test_shift_is_bits_per_symbol(alphabet_size, expected):
    assert misc.shift(alphabet_size) == expected


@pytest.mark.parametrize(
    "alphabet_size, expected", [(2, 1), (4, 3), (5, 7), (8, 7), (16, 15)]
)
def test_mask_covers_one_symbol(alphabet_size, expected):
    assert misc.mask(alphabet_size) == expected


@pytest.mark.parametrize("rows, expected", [(4, 1), (16, 2), (64, 3), (256, 4)])
def test_q_gram_size(rows, expected):
    assert misc.q_gram_size(rows, 4) == expected


def test_rc_tuple_complements_and_reverses():
    # AC -> GT: A=0, C=1 packs to 0b0001; G=2, T=3 packs to 0b1011.
    assert misc.rc_tuple(0b0001, 4, 2) == 0b1011
    assert misc.rc_tuple(0b1011, 4, 2) == 0b0001


def test_rc_tuple_is_an_involution():
    for code in range(16):
        assert misc.rc_tuple(misc.rc_tuple(code, 4, 2), 4, 2) == code


@pytest.mark.parametrize(
    "seq, expected",
    [
        ("ACGT", [0, 4]),
        ("NNACGTNNACGTN", [2, 6, 8, 12]),
        ("", []),
        ("NNN", []),
        ("ACGTN", [0, 4]),
        ("NACGT", [1, 5]),
        ("acgtNACGT", [0, 4, 5, 9]),
    ],
)
def test_preprocess_seq_finds_scannable_regions(seq, expected):
    assert misc.preprocess_seq(seq, 4, build_map()) == expected


def test_preprocess_seq_treats_every_non_base_as_a_break():
    assert misc.preprocess_seq("AC-GT", 4, build_map()) == [0, 2, 3, 5]
