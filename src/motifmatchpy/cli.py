"""Command line interface.

``motifmatchpy scan`` is the counterpart of the upstream ``moods-dna.py``
script; ``threshold`` and ``info`` expose the matrix utilities that are
otherwise only reachable from Python.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections.abc import Iterator, Sequence
from contextlib import contextmanager

from . import __version__, tools
from . import io as _io
from .api import Hit, Motif, MotifScanner, read_motifs

__all__ = ["main", "build_parser"]

PROG = "motifmatchpy"

_FORMAT_SEPARATORS = {"csv": ",", "tsv": "\t"}

SCAN_COLUMNS = (
    "sequence",
    "motif",
    "pos",
    "strand",
    "score",
    "match",
    "variant_match",
)

BED_COLUMNS = ("chrom", "start", "end", "name", "score", "strand")

INFO_COLUMNS = (
    "motif",
    "length",
    "order",
    "rows",
    "columns",
    "max_score",
    "min_score",
)


# ----------------------------------------------------------------------
# argument parsing
# ----------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=PROG,
        description="Match position weight matrices against DNA sequences.",
    )
    parser.add_argument("--version", action="version", version=f"{PROG} {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")

    # -- scan ----------------------------------------------------------
    scan_p = sub.add_parser(
        "scan",
        help="find motif occurrences in sequences",
        description="Find motif occurrences in DNA sequences.",
    )
    _add_matrix_arguments(scan_p)
    scan_p.add_argument(
        "-s",
        "--sequences",
        metavar="FILE",
        nargs="+",
        default=[],
        required=True,
        help="FASTA or plain-text sequence files (may be gzipped)",
    )

    cutoff = scan_p.add_argument_group("score cutoff (exactly one required)")
    cutoff.add_argument(
        "-p", "--p-value", metavar="P", type=float,
        help="report matches whose false-positive rate is at most P",
    )
    cutoff.add_argument(
        "-t", "--threshold", metavar="T", type=float,
        help="report matches scoring at least T",
    )
    cutoff.add_argument(
        "-B", "--best-hits", metavar="N", type=int,
        help="report approximately the N best-scoring matches per motif",
    )

    out = scan_p.add_argument_group("output")
    out.add_argument("-o", "--output", metavar="FILE", help="write here instead of stdout")
    out.add_argument(
        "-f", "--format", choices=("csv", "tsv", "bed"), default="csv",
        help="output format (default: csv)",
    )
    out.add_argument("--sep", metavar="S", help="override the field separator")
    out.add_argument("--header", action="store_true", help="write a header line")

    behaviour = scan_p.add_argument_group("search behaviour")
    behaviour.add_argument(
        "-R", "--no-rc", action="store_true",
        help="search the given strand only",
    )
    behaviour.add_argument(
        "--no-snps", action="store_true",
        help="ignore IUPAC codes standing for several nucleotides",
    )
    behaviour.add_argument(
        "--max-hits", metavar="N", type=int,
        help="give up on a motif once it exceeds N matches in a sequence",
    )
    behaviour.add_argument(
        "--batch", action="store_true",
        help="use --bg for every sequence instead of re-estimating the "
             "background, and so the thresholds, per sequence",
    )
    behaviour.add_argument(
        "--bg", metavar=("pA", "pC", "pG", "pT"), nargs=4, type=float,
        default=[0.25, 0.25, 0.25, 0.25],
        help="background for turning p-values into thresholds (default: uniform)",
    )
    behaviour.add_argument(
        "--threshold-precision", metavar="X", type=float,
        default=tools.DEFAULT_DP_PRECISION,
        help=f"precision of the p-value computation (default: {tools.DEFAULT_DP_PRECISION:g})",
    )
    scan_p.add_argument(
        "-v", "--verbose", action="count", default=0,
        help="report progress on stderr (-vv for more)",
    )
    scan_p.set_defaults(func=cmd_scan)

    # -- threshold -----------------------------------------------------
    th_p = sub.add_parser(
        "threshold",
        help="print the score threshold for a p-value",
        description="Print the score threshold each matrix needs to reach a p-value.",
    )
    _add_matrix_arguments(th_p)
    th_p.add_argument(
        "-p", "--p-value", metavar="P", type=float, required=True,
        help="p-value to convert into a score threshold",
    )
    th_p.add_argument(
        "--bg", metavar=("pA", "pC", "pG", "pT"), nargs=4, type=float,
        default=[0.25, 0.25, 0.25, 0.25],
        help="background distribution (default: uniform)",
    )
    th_p.add_argument(
        "--threshold-precision", metavar="X", type=float,
        default=tools.DEFAULT_DP_PRECISION,
        help="precision of the p-value computation",
    )
    th_p.add_argument("-o", "--output", metavar="FILE", help="write here instead of stdout")
    th_p.add_argument(
        "-f", "--format", choices=("csv", "tsv"), default="csv",
        help="output format (default: csv)",
    )
    th_p.add_argument("--sep", metavar="S", help="override the field separator")
    th_p.add_argument("--header", action="store_true", help="write a header line")
    th_p.set_defaults(func=cmd_threshold)

    # -- info ----------------------------------------------------------
    info_p = sub.add_parser(
        "info",
        help="summarise matrix files",
        description="Print the shape and score range of each matrix file.",
    )
    _add_matrix_arguments(info_p)
    info_p.add_argument("-o", "--output", metavar="FILE", help="write here instead of stdout")
    info_p.add_argument(
        "-f", "--format", choices=("csv", "tsv"), default="csv",
        help="output format (default: csv)",
    )
    info_p.add_argument("--sep", metavar="S", help="override the field separator")
    info_p.add_argument("--header", action="store_true", help="write a header line")
    info_p.set_defaults(func=cmd_info)

    return parser


def _add_matrix_arguments(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("matrices (at least one required)")
    group.add_argument(
        "-m", "--matrices", metavar="FILE", nargs="+", default=[],
        help="count matrices (.pfm) or dinucleotide models (.adm), "
             "converted to log-odds scores before matching",
    )
    group.add_argument(
        "-S", "--score-matrices", metavar="FILE", nargs="+", default=[],
        help="matrices that already hold scores, matched as they are",
    )
    group.add_argument(
        "--ps", metavar="P", type=float, default=0.01,
        help="pseudocount added per column in the log-odds conversion (default: 0.01)",
    )
    group.add_argument(
        "--lo-bg", metavar=("pA", "pC", "pG", "pT"), nargs=4, type=float,
        default=[0.25, 0.25, 0.25, 0.25],
        help="background used in the log-odds conversion (default: uniform)",
    )
    group.add_argument(
        "--log-base", metavar="X", type=float,
        help="logarithm base for the log-odds conversion (default: natural log)",
    )


# ----------------------------------------------------------------------
# shared helpers
# ----------------------------------------------------------------------
class UsageError(Exception):
    """A problem with the command line or the inputs it points at."""


def _load_motifs(args: argparse.Namespace) -> list[Motif]:
    if not args.matrices and not args.score_matrices:
        raise UsageError("no matrix files given (use -m or -S)")
    if args.log_base is not None and args.log_base <= 1:
        raise UsageError(f"--log-base must be greater than 1, got {args.log_base}")

    try:
        motifs = read_motifs(
            args.matrices,
            bg=args.lo_bg,
            pseudocount=args.ps,
            log_base=args.log_base,
            log_odds=True,
        )
        motifs += read_motifs(args.score_matrices, log_odds=False)
    except (OSError, ValueError, RuntimeError) as exc:
        raise UsageError(str(exc)) from exc
    return motifs


def _separator(args: argparse.Namespace) -> str:
    if getattr(args, "sep", None) is not None:
        return args.sep
    return _FORMAT_SEPARATORS.get(args.format, "\t")


@contextmanager
def _output(path: str | None) -> Iterator:
    if path is None:
        yield sys.stdout
        return
    try:
        handle = open(path, "w")  # noqa: SIM115 -- closed by this context manager
    except OSError as exc:
        raise UsageError(f"could not open {path} for writing: {exc}") from exc
    try:
        yield handle
    finally:
        handle.close()


def _log(args: argparse.Namespace, level: int, message: str) -> None:
    if getattr(args, "verbose", 0) >= level:
        print(f"{PROG}: {message}", file=sys.stderr)


# ----------------------------------------------------------------------
# commands
# ----------------------------------------------------------------------
def cmd_scan(args: argparse.Namespace) -> int:
    cutoffs = [args.p_value is not None, args.threshold is not None, args.best_hits is not None]
    if sum(cutoffs) == 0:
        raise UsageError("no score cutoff given (use -p, -t or -B)")
    if sum(cutoffs) > 1:
        raise UsageError("give only one of -p, -t and -B")
    if args.p_value is not None and not 0 < args.p_value <= 1:
        raise UsageError(f"-p must be in (0, 1], got {args.p_value}")
    if args.best_hits is not None and args.best_hits < 1:
        raise UsageError(f"-B must be at least 1, got {args.best_hits}")
    if args.batch and args.p_value is None:
        _log(args, 0, "warning: --batch only affects -p, ignoring it")

    motifs = _load_motifs(args)
    _log(args, 1, f"read {len(motifs)} matrices")

    # Without --batch a p-value threshold is recomputed from each sequence's own
    # base composition, which is what moods-dna does by default.
    background: Sequence[float] | str
    if args.p_value is not None and not args.batch:
        background = "auto"
    else:
        background = args.bg

    scanner = MotifScanner(
        motifs,
        p_value=args.p_value,
        threshold=args.threshold,
        bg=background,
        both_strands=not args.no_rc,
        threshold_precision=args.threshold_precision,
    )

    written = 0
    with _output(args.output) as out:
        writer = _writer(out, args)
        if args.header:
            writer.writerow(BED_COLUMNS if args.format == "bed" else SCAN_COLUMNS)

        for path in args.sequences:
            _log(args, 1, f"reading {path}")
            try:
                records = _io.read_sequences(path)
                for name, seq in records:
                    _log(args, 1, f"scanning {name} ({len(seq)} bp)")
                    if args.best_hits is not None:
                        hits = scanner.scan_best_hits(
                            seq, args.best_hits, sequence_name=name
                        )
                    else:
                        hits = scanner.scan(
                            seq,
                            sequence_name=name,
                            max_hits=args.max_hits,
                            include_variants=not args.no_snps,
                        )
                    _log(args, 2, f"{len(hits)} matches in {name}")
                    writer.writerows(_hit_row(hit, seq, args.format) for hit in hits)
                    written += len(hits)
            except OSError as exc:
                raise UsageError(f"could not read sequence file {path}: {exc}") from exc

    _log(args, 1, f"{written} matches in total")
    return 0


def cmd_threshold(args: argparse.Namespace) -> int:
    if not 0 < args.p_value <= 1:
        raise UsageError(f"-p must be in (0, 1], got {args.p_value}")
    motifs = _load_motifs(args)
    with _output(args.output) as out:
        writer = _writer(out, args)
        if args.header:
            writer.writerow(("motif", "threshold", "max_score", "min_score"))
        for motif in motifs:
            threshold = motif.threshold_from_p(
                args.p_value, args.bg, precision=args.threshold_precision
            )
            writer.writerow(
                (
                    motif.name,
                    _num(threshold),
                    _num(motif.max_score),
                    _num(motif.min_score),
                )
            )
    return 0


def cmd_info(args: argparse.Namespace) -> int:
    motifs = _load_motifs(args)
    with _output(args.output) as out:
        writer = _writer(out, args)
        if args.header:
            writer.writerow(INFO_COLUMNS)
        for motif in motifs:
            writer.writerow(
                (
                    motif.name,
                    str(motif.length),
                    str(motif.order),
                    str(len(motif.matrix)),
                    str(motif.width),
                    _num(motif.max_score),
                    _num(motif.min_score),
                )
            )
    return 0


# ----------------------------------------------------------------------
# formatting
# ----------------------------------------------------------------------
def _writer(handle, args: argparse.Namespace):
    """A csv writer, so a comma or tab inside a name cannot break a row."""
    return csv.writer(handle, delimiter=_separator(args), lineterminator="\n")


def _num(value: float) -> str:
    return f"{value:.6g}"


def _hit_row(hit: Hit, seq: str, fmt: str) -> tuple[str, ...]:
    if fmt == "bed":
        return (
            hit.sequence_name,
            str(hit.pos),
            str(hit.end),
            hit.name,
            _num(hit.score),
            hit.strand,
        )
    return (
        hit.sequence_name,
        hit.name,
        str(hit.pos),
        hit.strand,
        _num(hit.score),
        hit.matched_sequence(seq),
        hit.variant_sequence(seq) if hit.variants else "",
    )


# ----------------------------------------------------------------------
# entry point
# ----------------------------------------------------------------------
def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "func", None) is None:
        parser.print_help(sys.stderr)
        return 1
    try:
        return args.func(args)
    except UsageError as exc:
        print(f"{PROG}: error: {exc}", file=sys.stderr)
        return 2
    except BrokenPipeError:
        return 0
    except KeyboardInterrupt:
        print(f"{PROG}: interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
