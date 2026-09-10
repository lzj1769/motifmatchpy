"""The small value types the scanners return."""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ["Match", "Variant", "MatchWithVariant"]


@dataclass(frozen=True, slots=True)
class Match:
    """A motif occurrence: where it starts, and what it scored."""

    pos: int
    score: float

    def __iter__(self):
        yield self.pos
        yield self.score

    def __len__(self) -> int:
        return 2

    def __repr__(self) -> str:
        return f"Match(pos={self.pos}, score={self.score:.12g})"


@dataclass(frozen=True, slots=True, order=True)
class Variant:
    """A sequence edit: ``[start_pos, end_pos)`` replaced by ``modified_seq``.

    ``end_pos == start_pos + 1`` with a one-character replacement is a
    substitution, ``end_pos == start_pos`` an insertion, and an empty
    ``modified_seq`` a deletion.
    """

    start_pos: int
    end_pos: int
    modified_seq: str = ""

    def __repr__(self) -> str:
        return (
            f"Variant(start_pos={self.start_pos}, end_pos={self.end_pos}, "
            f"modified_seq={self.modified_seq!r})"
        )


@dataclass(frozen=True, slots=True)
class MatchWithVariant:
    """An occurrence that exists only once some variants are applied.

    ``pos`` is relative to the modified sequence, and ``variants`` holds indices
    into the variant list that was handed to the scanner.
    """

    pos: int
    score: float
    variants: tuple[int, ...] = field(default=())

    def __repr__(self) -> str:
        return (
            f"MatchWithVariant(pos={self.pos}, score={self.score:.12g}, "
            f"variants={list(self.variants)})"
        )
