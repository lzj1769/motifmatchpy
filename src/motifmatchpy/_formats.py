"""Content detection and DNA MEME/TRANSFAC count-matrix readers."""

from __future__ import annotations

import math
import re

from .io import open_text


def detect_format(filename) -> str:
    """Identify a motif format by its first meaningful line, not its suffix."""
    with open_text(filename) as handle:
        for line in handle:
            line = line.strip().lstrip("\ufeff")
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                return "jaspar"
            if line.startswith("MEME version"):
                return "meme"
            if line.split()[0] in {"AC", "ID", "VV", "XX", "DE", "P0", "PO", "//"}:
                return "transfac"
            try:
                float(line.split()[0])
            except ValueError as exc:
                raise ValueError(f"{filename}: unrecognized motif format") from exc
            return "numeric"
    raise ValueError(f"{filename}: empty motif file")


def _lines(filename):
    with open_text(filename) as handle:
        return [line.strip().lstrip("\ufeff") for line in handle]


def _append(result, name, positions, filename):
    if not name or not positions:
        raise ValueError(f"{filename}: missing motif name or matrix")
    if any(len(row) != 4 for row in positions):
        raise ValueError(f"{filename}: {name}: expected four DNA values per position")
    if any(not math.isfinite(v) or v < 0 for row in positions for v in row):
        raise ValueError(f"{filename}: {name}: values must be finite and non-negative")
    if any(existing == name for existing, _ in result):
        raise ValueError(f"{filename}: duplicate motif name {name!r}")
    result.append((name, [list(row) for row in zip(*positions, strict=True)]))


def meme(filename):
    """Read DNA MEME text probabilities, scaled by nsites (default 20)."""
    lines = _lines(filename)
    result = []
    name = None
    found = False
    i = 0
    while i < len(lines):
        line = lines[i]
        i += 1
        if line.startswith("ALPHABET"):
            if "=" not in line or "".join(line.split("=", 1)[1].split()) != "ACGT":
                raise ValueError(f"{filename}: only the MEME ACGT alphabet is supported")
        elif line.startswith("MOTIF ") or line.startswith("MOTIF\t"):
            if name is not None and not found:
                raise ValueError(f"{filename}: {name}: missing letter-probability matrix")
            name = "_".join(line.split()[1:])
            found = False
        elif line.startswith("letter-probability matrix:"):
            if name is None or found:
                raise ValueError(f"{filename}: unexpected letter-probability matrix")
            fields = dict(re.findall(r"(\w+)\s*=\s*(\S+)", line))
            try:
                width = int(fields["w"]) if "w" in fields else None
                alphabet = int(fields.get("alength", 4))
                sites = float(fields.get("nsites", 20))
            except ValueError as exc:
                raise ValueError(f"{filename}: {name}: invalid MEME matrix dimensions") from exc
            if alphabet != 4 or (width is not None and width < 1):
                raise ValueError(f"{filename}: {name}: expected a positive-width DNA matrix")
            if not math.isfinite(sites) or sites <= 0:
                raise ValueError(f"{filename}: {name}: nsites must be positive and finite")
            positions = []
            while i < len(lines):
                row = lines[i].split()
                if not row:
                    break
                try:
                    values = [float(v) for v in row]
                except ValueError:
                    break
                if len(values) != 4 or any(not math.isfinite(v) or not 0 <= v <= 1
                                           for v in values):
                    raise ValueError(f"{filename}: {name}: invalid probability row")
                if not math.isclose(sum(values), 1.0, abs_tol=1e-4):
                    raise ValueError(f"{filename}: {name}: probability row must sum to 1")
                positions.append([v * sites for v in values])
                i += 1
            if width is not None and len(positions) != width:
                raise ValueError(f"{filename}: {name}: expected {width} probability rows")
            _append(result, name, positions, filename)
            found = True
    if name is not None and not found:
        raise ValueError(f"{filename}: {name}: missing letter-probability matrix")
    if not result:
        raise ValueError(f"{filename}: no MEME probability matrices found")
    return result


def transfac(filename):
    """Read TRANSFAC DNA counts; transpose position rows into A/C/G/T rows."""
    result = []
    accession = identifier = None
    columns = None
    positions = []
    for line in _lines(filename):
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        key = fields[0]
        if key == "//":
            if accession is not None or identifier is not None or columns is not None:
                parts = [p for p in (accession, identifier) if p]
                name = "_".join(dict.fromkeys(parts))
                _append(result, name, positions, filename)
            accession = identifier = columns = None
            positions = []
        elif key in {"AC", "ID"}:
            value = "_".join(fields[1:]).rstrip(";")
            if not value or (key == "AC" and accession is not None) or (
                key == "ID" and identifier is not None
            ):
                raise ValueError(f"{filename}: invalid or repeated {key}; records need //")
            if key == "AC":
                accession = value
            else:
                identifier = value
        elif key in {"P0", "PO"}:
            if columns is not None or len(fields) != 5 or set(fields[1:]) != set("ACGT"):
                raise ValueError(f"{filename}: expected a TRANSFAC P0 A C G T header")
            columns = fields[1:]
        elif key.isdigit():
            if columns is None or int(key) != len(positions) + 1 or len(fields) not in {5, 6}:
                raise ValueError(f"{filename}: invalid TRANSFAC position row: {line}")
            try:
                values = [float(v) for v in fields[1:5]]
            except ValueError as exc:
                raise ValueError(f"{filename}: invalid TRANSFAC counts") from exc
            positions.append([values[columns.index(base)] for base in "ACGT"])
    if accession is not None or identifier is not None or columns is not None:
        raise ValueError(f"{filename}: unterminated TRANSFAC record (expected //)")
    if not result:
        raise ValueError(f"{filename}: no TRANSFAC motifs found")
    return result
