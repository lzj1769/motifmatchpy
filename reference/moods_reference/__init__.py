"""The original MOODS C++ core, exposed for comparison in the test suite.

``motifmatchpy`` is a pure-Python reimplementation of MOODS. This package wraps
the C++ sources it was ported from, so the tests can assert that the two agree
on real data rather than only on hand-computed examples.

The vendored C++ carries a small number of fixes without which it cannot serve
as an oracle -- it would crash or silently drop matches. Each one is recorded in
``core/VENDORED.md`` together with the upstream behaviour it replaces, and
``tests/test_reference_parity.py`` pins the differences that remain observable.
"""

from ._core import (  # noqa: F401
    MOODS_VERSION,
    Match,
    MatchWithVariant,
    Variant,
    misc,
    parsers,
    scan,
    tools,
)

__all__ = [
    "MOODS_VERSION",
    "Match",
    "MatchWithVariant",
    "Variant",
    "misc",
    "parsers",
    "scan",
    "tools",
]
