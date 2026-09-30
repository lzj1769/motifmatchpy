# threshold

Calculate each matrix's score cutoff under a background model without scanning
a sequence. `-p/--p-value` is required and must be in `(0, 1]`.

```sh
motifmatchpy threshold -S docs/examples/toy.scores -p 0.02 --header
```

```text
motif,threshold,max_score,min_score
toy.scores,3.0005,6,-3
```

## Options

Use the [shared matrix/output options](index.md), plus:

| Option | Meaning | Default |
| --- | --- | --- |
| `-p`, `--p-value P` | Tail-probability cutoff under the background model | Required |
| `--bg pA pC pG pT` | Threshold background in A/C/G/T order | Uniform |
| `--threshold-precision X` | Positive score discretization precision | `2000` |

Larger precision increases the dynamic-programming state space and may change
the approximate cutoff, using more time and memory. A p-value is not an adjusted
false discovery rate, and a computed threshold can exclude all positions.

## Count matrix and custom background

```sh
motifmatchpy threshold -m docs/examples/toy.pfm -p 0.02 \
  --lo-bg 0.3 0.2 0.2 0.3 --bg 0.3 0.2 0.2 0.3 \
  --threshold-precision 10000 -f tsv --header -o thresholds.tsv
```

The columns are `motif`, `threshold`, `max_score`, `min_score`. Scores and
thresholds use the input matrix's units, or the logarithm base used to convert
counts. To reproduce this cutoff during a scan, use the same matrix conversion,
`-p`, `--bg`, and precision settings **with `--batch`**.

The command reports the supplied motif's threshold only. A reverse-complement
motif may have a different cutoff under an asymmetric background; the scanner
calculates its strand-specific thresholds internally.
