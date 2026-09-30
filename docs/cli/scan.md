# scan

Find motif occurrences in one or more FASTA or plain-text files, including gzip
input. Select **exactly one** cutoff: `-p`, `-t`, or `-B`.
Alternatively, provide `--regions BED --genome FASTA` to scan genomic intervals;
this defaults to BED output and `-p 1e-4`. See the
[ATAC-seq peak workflow](../peaks.md) for a full JASPAR example.

```sh
motifmatchpy scan -S docs/examples/toy.scores -s docs/examples/toy.fa \
  -t 6 --header
```

This returns positions 0 and 6 on `+`, and position 1 on `-`.
See the [quick start](../quickstart.md) for the exact CSV output.

## Options

In addition to the [shared matrix/output options](index.md):

| Option | Behavior |
| --- | --- |
| `-s`, `--sequences FILE ...` | Sequence files; mutually exclusive with peak mode |
| `--regions BED` / `--genome FASTA` | Peak intervals and reference genome; provide both |
| `-p`, `--p-value P` | Derive a score cutoff; `0 < P <= 1` |
| `-t`, `--threshold T` | Report scores at least T, using the matrix's score units |
| `-B`, `--best-hits N` | Search for approximately N matches per motif and strand, per record |
| `-R`, `--no-rc` | Search only the supplied strand |
| `--no-snps` | Disable extra matches from IUPAC alternatives |
| `--max-hits N` | Stop collecting ordinary hits after N per motif/strand; keep the first hits |
| `--batch` | Keep `--bg` fixed for p-value thresholds across records |
| `--bg pA pC pG pT` | Threshold background; uniform by default |
| `--threshold-precision X` | Dynamic-programming discretization precision; default 2000 |
| `-v`, `--verbose` | Write progress to stderr; `-vv` includes match counts |

The default search checks both strands. With `-p` or `-t`, IUPAC alternatives
are included unless `--no-snps` is set. Unknown bases break ordinary matches.
`N` is not expanded into all four bases; see [file formats](../formats.md).

## P-value scan

```sh
motifmatchpy scan -m docs/examples/toy.pfm -s docs/examples/toy.fa \
  -p 0.02 --batch --no-snps --header -o hits.csv
```

This converts counts to natural-log scores, uses a uniform background for both
conversion and threshold calculation, and writes three ordinary matches.
In `-s` mode without `--batch`, sequence composition changes the threshold for
each record. Peak mode always uses the fixed `--bg`.

## BED output and forward strand only

```sh
motifmatchpy scan -S docs/examples/toy.scores -s docs/examples/toy.fa \
  -t 6 -R -f bed -o hits.bed
```

```text
toy	0	3	toy.scores	6	+
toy	6	9	toy.scores	6	+
```

BED coordinates are zero-based and half-open. The fifth column contains the raw
motif score, **not** a normalized integer from 0 to 1000. Some BED consumers
require rescaling. Omit `--header` for conventional BED input to other tools.

## Approximate best hits

```sh
motifmatchpy scan -S docs/examples/toy.scores -s docs/examples/toy.fa \
  -B 1 --header
```

This searches a score threshold rather than sorting and truncating to exactly N
rows. Ties can yield extra hits, and very common motifs may be abandoned.
The target is applied independently to each searched strand, not to a pooled
top-N list. Best-hits estimates its own sequence background.

`-B` does not run the IUPAC variant scan and does not apply `--max-hits`.
With `-p` or `-t`, `--max-hits` limits ordinary hits, not additional variant hits.
For variant analysis, use `-p` or `-t`, or the low-level
[variant API](../advanced.md#explicit-variants).

## Output columns

In peak mode, positions include the peak start offset and refer to the reference
genome. In `-s` mode they are relative to the individual sequence record.

CSV/TSV columns, in order: `sequence`, `motif`, `pos`, `strand`, `score`,
`match`, `variant_match`. The last field is empty for ordinary hits. The matched
bases are shown in input-sequence orientation even for `-` hits.

BED columns: `chrom`, `start`, `end`, `name`, `score`, `strand`. BED omits variant
details; use CSV/TSV to retain them. Coordinates and interpretation are described
in [file formats](../formats.md).
