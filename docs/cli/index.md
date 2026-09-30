# Command line

The three subcommands are `scan`, `threshold`, and `info`. Each accepts one or
more matrix files. Commands in this guide run from the repository root.

```sh
motifmatchpy --version
motifmatchpy --help
motifmatchpy scan --help
motifmatchpy threshold --help
motifmatchpy info --help
```

## Matrix arguments shared by all subcommands

| Option | Meaning | Default |
| --- | --- | --- |
| `-m`, `--matrices FILE ...` | PFM counts or ADM terms; convert to log-odds | None |
| `-S`, `--score-matrices FILE ...` | Numeric score matrices; do not convert | None |
| `--ps P` | Pseudocount used during conversion | `0.01` |
| `--lo-bg pA pC pG pT` | Background for log-odds conversion | `0.25 0.25 0.25 0.25` |
| `--log-base X` | Logarithm base, greater than 1 | Natural logarithm |

Supply at least one file via `-m` or `-S`; both may be used together. Each file
contains one motif. Shell globs such as `-m matrices/*.pfm` are supported through
shell expansion. Motif names default to filenames including their suffixes.
Conversion options do not transform matrices passed through `-S`.

## Output options

All subcommands support `-o/--output FILE`, `-f/--format csv|tsv`, `--header`,
and `--sep S`. Output goes to stdout as CSV without a header unless specified.
`scan` additionally supports `-f bed`. `--sep` accepts a single separator
character and overrides the selected format's separator.

```sh
motifmatchpy info -m docs/examples/toy.pfm -f tsv --header -o motifs.tsv
```

Scores are printed with **six significant digits**. Use the Python API when
retaining full floating-point precision matters. Files passed to `-o` are
overwritten; parent directories must already exist.

## Two different backgrounds

- `--lo-bg` determines how input counts become log-odds scores.
- `--bg` determines the score distribution used for p-value thresholds.

For `scan -p`, the default is to estimate the threshold background separately
for each sequence record. Add `--batch` to use `--bg` for every record. Setting
`--bg` alone does **not** override this per-record estimation. `threshold` always
uses its `--bg` argument and has no `--batch` option.

```sh
motifmatchpy scan -m docs/examples/toy.pfm -s docs/examples/toy.fa \
  -p 0.02 --batch --bg 0.3 0.2 0.2 0.3 --lo-bg 0.3 0.2 0.2 0.3
```

Backgrounds follow A, C, G, T order; use positive probabilities summing to one.
For reproducible thresholds, record both backgrounds, pseudocount, logarithm
base, p-value, and threshold precision.

## Exit status

Success returns 0; handled argument or input errors return 2; an interrupt
returns 130. Calling `motifmatchpy` without a subcommand prints help and returns
1. A successful scan can have no matches. Use `--header` to distinguish an empty
result table from an absent output file.
