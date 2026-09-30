# motifmatchpy

**Find DNA motif occurrences from Python or the command line.** motifmatchpy
implements the MOODS motif-matching algorithms in Python, with Numba compiling
the scanning loops at runtime. It supports multiple motifs, both strands,
p-value thresholds, higher-order models, and sequence variants.

## Start here

- [Install](installation.md) with Python 3.12 or newer.
- [Run the quick start](quickstart.md) using the bundled tiny example files.
- Use [scan](cli/scan.md) to find matches, [threshold](cli/threshold.md) to convert
  p-values to score cutoffs, or [info](cli/info.md) to inspect matrices.
- Scan [specified genomic regions with a reference genome and motif files](peaks.md).
- Integrate a reusable scanner with the [Python guide](python.md) and
  [API reference](reference/api.md).

## Choose an interface

| Task | Interface |
| --- | --- |
| Scan FASTA and export CSV, TSV, or BED | `motifmatchpy scan` |
| Inspect motif dimensions and score ranges | `motifmatchpy info` |
| Calculate thresholds under a background distribution | `motifmatchpy threshold` |
| Track motif names, strands, and coordinates in Python | `Motif`, `MotifScanner`, `Hit` |
| Work with matrices, custom alphabets, or explicit indels | `tools`, `scan.Scanner`, `parsers`, `misc` |

Motif matches are sequence-level candidates, not evidence of binding in a cell.
The p-value cutoff refers to a background sequence model, not a multiple-testing
adjusted significance level.

## Version and provenance

These pages describe the repository's `main` branch. Check `motifmatchpy --version`
for your installed release; fixes on `main` can precede the next package release.
The algorithm baseline is MOODS 1.9.4.1. Known corrections to upstream behavior
are explained in [compatibility](compatibility.md).

The project includes the upstream GPLv3 and Biopython license texts in the
[repository](https://github.com/lzj1769/motifmatchpy). When describing the
underlying algorithms, cite Korhonen et al.,
[Fast motif matching revisited: high-order PWMs, SNPs and indels](https://doi.org/10.1093/bioinformatics/btw683),
*Bioinformatics* (2017). A separate motifmatchpy paper is not yet available;
record the package version and repository URL in reproducible analyses.
