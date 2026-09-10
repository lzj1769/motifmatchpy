"""Shared fixtures.

The example data ships with MOODS: 16 JASPAR count matrices, one adjacent
dinucleotide model, and 50 kb of human chromosome 1 (which contains long runs
of ``N``, so the region-splitting code gets exercised on real input).
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"
MATRICES = DATA / "matrices"
SEQUENCE_FILE = DATA / "seq" / "chr1-5k-55k.fa"


def pytest_collection_modifyitems(config, items):
    """Skip the C++ comparison tests when the reference build is not installed."""
    try:
        import moods_reference  # noqa: F401
    except ImportError:
        skip = pytest.mark.skip(
            reason="moods-reference is not installed; run `uv sync` to build it"
        )
        for item in items:
            if "reference" in item.keywords:
                item.add_marker(skip)


@pytest.fixture(scope="session")
def matrix_dir() -> Path:
    return MATRICES


@pytest.fixture(scope="session")
def pfm_files() -> list[Path]:
    return sorted(MATRICES.glob("*.pfm"))


@pytest.fixture(scope="session")
def pfm_file() -> Path:
    return MATRICES / "MA0001.pfm"


@pytest.fixture(scope="session")
def adm_file() -> Path:
    return MATRICES / "E2F3.adm"


@pytest.fixture(scope="session")
def sequence_file() -> Path:
    return SEQUENCE_FILE


@pytest.fixture(scope="session")
def dna() -> str:
    """The example chromosome fragment, mixed case, ``N`` runs and all."""
    lines = SEQUENCE_FILE.read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


@pytest.fixture(scope="session")
def sample(dna: str) -> str:
    """A 20 kb slice, small enough to brute-force scan in Python."""
    return dna[:20000]


@pytest.fixture
def rng() -> random.Random:
    return random.Random(20240611)


@pytest.fixture(scope="session")
def flat4() -> list[float]:
    return [0.25, 0.25, 0.25, 0.25]


@pytest.fixture(scope="session")
def skewed4() -> list[float]:
    """A deliberately lopsided background, where symmetry cannot hide a bug."""
    return [0.31, 0.19, 0.23, 0.27]


@pytest.fixture(scope="session")
def polya() -> list[list[float]]:
    """A toy matrix scoring 1 per A and 0 otherwise, so hits are obvious by eye."""
    return [[1.0] * 5, [0.0] * 5, [0.0] * 5, [0.0] * 5]
