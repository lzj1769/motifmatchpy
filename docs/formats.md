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

This is a **numeric table**, not every format named JASPAR/PFM. Remove row labels
such as `A [ ... ]`, brackets, and multi-motif database wrappers first. The parser
stops each line at the first nonnumeric token and drops lines with no numbers;
labels at the start of a row therefore discard that row.

## ADM models

A DNA `.adm` file has 20 nonempty numeric rows:

- 16 equal-length rows for adjacent pairs, ordered `AA, AC, AG, AT, CA, ... , TT`.
- Four rows for the initial A/C/G/T marginal terms; only the first value of each
  is used during log-odds conversion.

An ADM with W conditional columns spans W+1 sequence positions. The reader
chooses a parser from the suffix first, then tries the other supported format.
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
