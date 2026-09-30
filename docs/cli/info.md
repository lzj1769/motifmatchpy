# info

Inspect motif dimensions, model order, and achievable score range.
No sequence or score cutoff is required.

```sh
motifmatchpy info -S docs/examples/toy.scores --header
```

```text
motif,length,order,rows,columns,max_score,min_score
toy.scores,3,0,4,3,6,-3
```

## Inspect a converted count matrix

```sh
motifmatchpy info -m docs/examples/toy.pfm --log-base 2 --header -f tsv
```

`-m` computes log-odds scores before reporting the range; `--log-base 2` gives
scores in bits. Use `-S` when the file already contains scores. All
[shared matrix/output options](index.md) apply; this subcommand has no additional
search options.

## Interpret the columns

| Column | Meaning |
| --- | --- |
| `motif` | Name, normally the input filename |
| `length` | Number of sequence bases spanned |
| `order` | 0 for a PWM; 1 for adjacent dinucleotides; higher for longer contexts |
| `rows` | Number of scoring rows; `4 ** (order + 1)` for DNA |
| `columns` | Number of matrix columns |
| `max_score` | Highest achievable score |
| `min_score` | Lowest achievable score |

For a higher-order model, `length = columns + order`. Overlapping contexts
constrain the achievable score range; it is not generally the independent sum
of column maxima/minima.
