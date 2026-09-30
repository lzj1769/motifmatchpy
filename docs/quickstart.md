# Quick start

Run these commands from the repository root after [installing](installation.md).
Alternatively, download [toy.fa](examples/toy.fa), [toy.pfm](examples/toy.pfm),
and [toy.scores](examples/toy.scores) and adjust the paths below.

`toy.fa` contains `ACGTTTACG`. The three-position motif favors `ACG`.
The score matrix awards 2 for each preferred base and -1 otherwise, so a perfect
match scores 6. Using explicit scores makes the example independent of p-value
discretization and background estimation.

## 1. Inspect a motif

```sh
motifmatchpy info -S docs/examples/toy.scores --header
```

```text
motif,length,order,rows,columns,max_score,min_score
toy.scores,3,0,4,3,6,-3
```

## 2. Calculate a threshold

```sh
motifmatchpy threshold -S docs/examples/toy.scores -p 0.02 --header
```

The default uniform background assigns probability `1/64` to a particular
three-base word. The discretized cutoff lies just above 3, so only perfect
matches scoring 6 qualify for this matrix.

```text
motif,threshold,max_score,min_score
toy.scores,3.0005,6,-3
```

## 3. Scan both strands

```sh
motifmatchpy scan -S docs/examples/toy.scores -s docs/examples/toy.fa -t 6 --header
```

```text
sequence,motif,pos,strand,score,match,variant_match
toy,toy.scores,0,+,6,ACG,
toy,toy.scores,1,-,6,CGT,
toy,toy.scores,6,+,6,ACG,
```

Positions are zero-based. The reverse-strand match at position 1 covers `CGT`
in the supplied sequence; its reverse complement is `ACG`.

## 4. Use count matrices

```sh
motifmatchpy scan -m docs/examples/toy.pfm -s docs/examples/toy.fa \
  -p 0.02 --batch --no-snps --header
```

`-m` converts counts into log-odds scores; `-S` uses scores directly.
`--batch` keeps the uniform background fixed for threshold calculation.
See the [CLI guide](cli/index.md) for the distinct conversion and threshold
backgrounds, and the [Python guide](python.md) for the same workflow in code.
