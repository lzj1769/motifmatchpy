# motifmatchpy

Position weight matrix matching for DNA: a pure-Python port of
[MOODS](https://github.com/jhkorhonen/MOODS), with the inner loops compiled by
[numba](https://numba.pydata.org/).

Scans DNA for occurrences of position weight matrices — hundreds of matrices
against megabases of sequence in well under a second — and handles first- and
higher-order models, custom alphabets, p-value thresholds, and the matches that
only appear once sequence variants are applied.

Results are **identical to the original C++ implementation, bit for bit**. The
test suite checks that against a build of it kept in `reference/`.

## Install

```sh
uv add motifmatchpy      # or: pip install motifmatchpy
```

Python 3.12 or newer. No compiler needed; the only dependencies are numpy and
numba.

## Python API

```python
import motifmatchpy as mm

motifs = mm.read_motifs(["MA0001.pfm", "MA0002.pfm", "E2F3.adm"])
scanner = mm.MotifScanner(motifs, p_value=1e-4)

for hit in scanner.scan_fasta("chr1.fa"):
    print(hit.sequence_name, hit.name, hit.pos, hit.strand, round(hit.score, 3))
```

A `Motif` knows its name, length and model order; a `Hit` knows where it is,
which strand it is on, and what it scored:

```python
motif = mm.Motif.from_file("MA0001.pfm")
motif.length, motif.order, motif.max_score
motif.threshold_from_p(1e-4)
motif.reverse_complement()

hits = scanner.scan(sequence, sequence_name="chr1")
hits[0].matched_sequence(sequence)      # the bases it covers
scanner.count(sequence)                 # occurrences per motif, both strands
```

`MotifScanner` takes exactly one of `p_value` or `threshold` (or neither, if you
only want `scan_best_hits`), scans both strands unless told otherwise, and
re-estimates the background per sequence with `bg="auto"`:

```python
mm.MotifScanner(motifs, p_value=1e-4, bg="auto")        # threshold per sequence
mm.MotifScanner(motifs, threshold=8.0, both_strands=False)
mm.MotifScanner(motifs).scan_best_hits(sequence, 1000)  # no threshold needed
```

### Sequence variants

Ambiguity codes and indels can create matches the plain sequence does not have.
`include_variants` reports them, each carrying the variants responsible:

```python
hits = scanner.scan(sequence, include_variants=True)
for hit in hits:
    if hit.variants:
        print(hit.pos, hit.variant_sequence(sequence), hit.variants)
```

### Working in bare matrices

`motifmatchpy.tools`, `.scan`, `.parsers` and `.misc` mirror the MOODS modules
of the same names, so code written against MOODS ports by changing the import:

```python
from motifmatchpy import parsers, scan, tools

bg = tools.bg_from_sequence_dna(seq, 1)
matrix = parsers.pfm_to_log_odds("MA0001.pfm", bg, 0.01)
threshold = tools.threshold_from_p(matrix, bg, 1e-4)
matches = scan.scan_dna(seq, [matrix], bg, [threshold])[0]
```

## Command line

```sh
motifmatchpy scan -m matrices/*.pfm -s chr1.fa -p 1e-5
motifmatchpy scan -m matrices/*.pfm -s chr1.fa -p 1e-5 -f bed -o hits.bed
motifmatchpy scan -m MA0001.pfm -s chr1.fa -B 1000        # ~1000 best matches
motifmatchpy threshold -m matrices/*.pfm -p 1e-4
motifmatchpy info -m matrices/*.pfm
```

`scan` writes CSV by default (`-f tsv`, `-f bed` also available), reads FASTA or
plain text, gzipped or not, and takes exactly one of `-p` (p-value), `-t`
(absolute score) or `-B` (approximate number of best matches). `-R` searches the
given strand only, `--no-snps` ignores ambiguity codes, and `--batch` keeps the
thresholds fixed instead of re-deriving them from each sequence. `motifmatchpy
scan -h` lists the rest.

## File formats

`.pfm`
: A JASPAR-style count matrix: one whitespace-separated row per base (A, C, G,
  T), one column per motif position.

`.adm`
: An adjacent dinucleotide model: 16 rows of first-order conditional terms
  followed by 4 rows of zero-order terms for the first position.

The format is taken from the suffix where possible and otherwise detected by
trying both, so a mislabelled file still reads.

## Performance

2 Mb of human chromosome 1, 32 matrices (16 JASPAR profiles and their reverse
complements), on an M-series Mac:

| threshold | matches | motifmatchpy | MOODS C++ |
|---|---:|---:|---:|
| p = 1e-4 | 12,760 | 0.040 s | 0.014 s |
| p = 1e-3 | 88,184 | 0.117 s | 0.036 s |
| p = 1e-2 | 677,356 | 0.509 s | 0.160 s |
| p = 1e-1 | 5,618,812 | 2.837 s | 0.835 s |

Roughly three times the C++ runtime, for the same matches. The first scan in a
process also pays a one-off second or so to JIT-compile the kernels; numba
caches the result on disk, so later runs skip it.

Reuse a `MotifScanner` or a `scan.Scanner` across sequences — preparing the
matrices costs more than scanning a short sequence does.

## How it works

Scanning position by position against every matrix would be far too slow, so
MOODS — and this port — read a fixed-width *window* of the sequence, pack its
contents into an integer, and look that integer up in a table listing the motifs
that could still reach their threshold from there. Motifs no longer than the
window are settled by the lookup alone. Longer ones are then scored outwards
from the window, one position at a time, in order of how much each position
discriminates, and abandoned as soon as the best possible completion falls
short.

`src/motifmatchpy/_motif.py` builds those tables, `_kernels.py` holds the
compiled scanning loop, and `scan.py` ties them together.

## Development

```sh
uv sync                          # installs the package and the C++ reference
uv run pytest                    # 521 tests, including the C++ parity suite
uv run pytest -m "not slow"      # skip the exhaustive sweeps
uv run pytest -m reference       # only the comparisons against MOODS
```

`reference/` holds the original MOODS C++ core with nanobind bindings, built as
a separate development-only package. It exists so the tests can assert that this
implementation and the original agree exactly, on real data and on tens of
thousands of generated cases.

Porting the algorithms turned up several genuine bugs in the C++ original —
including matches silently dropped at the edge of every `N` run, and wrong
p-value thresholds under a non-uniform background. `motifmatchpy` implements the
corrected behaviour; each divergence is documented, and patched in the reference
copy so the parity tests still compare like with like. See
[`reference/core/VENDORED.md`](reference/core/VENDORED.md).

## Licence

GPL v3 or later, or the Biopython licence — the same dual licence as MOODS,
whose algorithms this implements. See `COPYING.GPLv3` and `COPYING.BIOPYTHON`.

MOODS is by Pasi Rastas, Janne H. Korhonen and Petri Martinmäki:
[Fast motif matching revisited: high-order PWMs, SNPs and indels](https://doi.org/10.1093/bioinformatics/btw683),
*Bioinformatics* 2017.
