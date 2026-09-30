# File formats and coordinates

## PFM counts and score matrices

Use four whitespace-separated numeric rows in A, C, G, T order, one column per
motif position. Each row must have the same positive number of columns.
The [example PFM](examples/toy.pfm) is:

```text
10 0 0
0 10 0
0 0 10
0 0 0
```

`-m` or `Motif.from_file()` interprets these as counts/frequencies. Conversion
adds a background-weighted pseudocount, normalizes each column, divides by the
background probability, and takes a logarithm. Natural logarithms are default.
`-S` or `log_odds=False` interprets numeric entries as scores directly.

`-m/--matrices` and Python's `read_motifs()` detect formats by content:
`>` headers identify JASPAR, `MEME version` identifies MEME text, TRANSFAC uses
record tags, and numeric tables are distinguished by their shape. File suffixes
need not match. All formats support gzip. Once a collection format is recognized,
malformed records raise errors rather than falling back to numeric-table parsing.

Use `read_motifs()` for collections. `Motif.from_file()` accepts one motif and
rejects multi-motif input instead of silently selecting the first record.

## JASPAR collections

Use `-m FILE` for JASPAR download files with one or many motifs:

```text
>MA0000.1 TEST
A [10 0 0]
C [0 10 0]
G [0 0 10]
T [0 0 0]
```

Each record has a `>ID TF_NAME` header and four rows in A/C/G/T order.
Unlabelled numeric rows following a header (JASPAR's PFM export) are also
supported. Counts must be finite, non-negative, and equal-length. Empty lines
and `#` comments are ignored; gzip is supported. Duplicate combined ID/name
headers in one collection are rejected. `-m` always converts counts;
use `-S` for precomputed numeric scores.

For genomic BED3/narrowPeak inputs and output coordinates, see the
[genomic region workflow](peaks.md).

## MEME text collections

DNA MEME minimal text files start with `MEME version` and contain `MOTIF ID NAME`
records with `letter-probability matrix:` blocks. Matrix rows are positions,
columns are A/C/G/T. Probabilities must be finite, in [0, 1], and sum to one
within rounding tolerance. The parser checks `w` and `alength` when supplied;
without `w`, it infers the number of consecutive probability rows.

Probabilities are multiplied by `nsites` to obtain effective counts before the
usual pseudocount/log-odds conversion. Missing `nsites` defaults to 20, following
[MEME's format specification](https://meme-suite.org/meme/doc/meme-format.html).
For example, `nsites=10` and a probability of 0.7 produce a count of 7.
File background and strand declarations do not override `--lo-bg`, `--bg`, or
`-R`. Protein/custom alphabets, MEME XML/HTML, and files containing only log-odds
matrices are not supported by this count-input reader.

## TRANSFAC collections

TRANSFAC records contain `AC` and/or `ID` identifiers, a `P0 A C G T` (or `PO`)
header, consecutively numbered position rows starting at 1, and a terminating
`//`. Each row contains four finite, non-negative counts and optionally a
consensus symbol. Column permutations in the header are reordered into A/C/G/T.
Other metadata fields are ignored. Names combine AC and ID when both exist.
The [TRANSFAC format example](https://biopython.org/docs/latest/Tutorial/chapter_motifs.html)
illustrates the record layout.

## ADM models

A DNA `.adm` file has 20 nonempty numeric rows:

- 16 equal-length rows for adjacent pairs, ordered `AA, AC, AG, AT, CA, ... , TT`.
- Four rows for the initial A/C/G/T marginal terms; only the first value of each
  is used during log-odds conversion.

An ADM with W conditional columns spans W+1 sequence positions. Numeric table shape determines whether the data is PFM or ADM.
For raw ADM scores (`-S`), the 16 conditional rows are used directly; the four
trailing rows are required by the file parser but are not folded into scores.

File loading supports zero-order PFM and first-order ADM models. Construct
second- and higher-order `Motif` matrices directly in Python; see
[advanced models](advanced.md).

## Sequence files

FASTA records start with `>`. By default, the record name is the header's first
whitespace-delimited token. Python's `read_fasta` and `read_sequences` accept
`full_header=True` to retain the entire header. A plain-text file represents one
sequence and is named after its filename. Gzip is detected by file contents.

Input is ASCII; uppercase and lowercase DNA are equivalent. Non-ACGT characters
break ordinary matches. IUPAC alternatives `W R S Y K M B D H V` can generate
variant matches. `N` and unrecognized characters remain unknown and are not
expanded into alternatives. Whitespace at line ends is stripped; avoid spaces
inside sequence lines.

## Ordinary coordinates and strands

All positions are zero-based. A motif of length L at `pos` covers
`[pos, pos + L)`, regardless of strand. Reverse hits are reported relative to
the original input sequence; the scanner reverses the motif, not the sequence.

For `ACGTTTACG`, the reverse-strand hit has `pos=1`, `end=4`, and input slice
`CGT`. CSV's `match` and `Hit.matched_sequence()` retain this input orientation.
BED uses the same half-open interval.

## Variant coordinates

`Variant(start_pos, end_pos, modified_seq)` replaces the original interval
`[start_pos, end_pos)`. Equal start/end indicates an insertion; an empty
replacement indicates a deletion. `MatchWithVariant.pos` is relative to the
modified sequence. Its `variants` tuple contains indices into the supplied
variant list. This is not a VCF reader or a reference-coordinate liftover tool.

For substitutions the coordinate system is unchanged. For indels, retain the
edits alongside the result and map coordinates explicitly before intersecting
with reference-genome annotations.
