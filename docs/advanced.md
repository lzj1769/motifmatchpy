# Advanced models

The low-level API operates on score matrices and returns one result group per
matrix. It does not automatically add reverse complements or motif names.

## Bare-matrix scanning

```python
from motifmatchpy import parsers, scan, tools

bg = tools.flat_bg(4)
matrix = parsers.pfm("docs/examples/toy.scores")
threshold = tools.threshold_from_p(matrix, bg, 0.02)
scanner = scan.Scanner(window_size=7)
scanner.set_motifs([matrix], bg, [threshold])
assert [(m.pos, m.score) for m in scanner.scan("ACGTTTACG")[0]] == [(0, 6), (6, 6)]
```

Reuse `Scanner` for multiple sequences. `scan.scan_dna()` creates a scanner for
a one-off call. `scan.naive_scan_dna()` checks positions directly and is useful
for small independent comparisons. For a higher-order matrix, pass `a=4` to
utility functions such as `threshold_from_p`, `max_score`, `reverse_complement`,
and `naive_scan_dna` so they interpret overlapping contexts correctly.

## Higher-order count conversion

A DNA order-k matrix has `4 ** (k + 1)` rows. Its columns score overlapping
(k+1)-base words in lexicographic A/C/G/T order. A W-column matrix spans W+k
bases. Lower-order terms supply the leading positions that lack full context.

```python
import motifmatchpy as mm

# A simple uniform second-order model: 64 triplet rows, 3 columns.
conditional = [[1.0, 1.0, 1.0] for _ in range(64)]
low = [[1.0] * 4, [1.0] * 16]
scores = mm.tools.log_odds_high_order(conditional, low, [0.25] * 4, 0.01, a=4)
motif = mm.Motif("uniform-order2", scores)
assert (motif.order, motif.width, motif.length) == (2, 3, 5)
```

Lower-order row r requires at least `4 ** (r + 1)` values; rows need not have
equal lengths. A zero-order model uses `low_order_terms=[]`. Pass `log_base=2`
for bits. Use `log_odds_high_order` explicitly when requesting a log base;
the overloaded `log_odds` spelling does not expose every C++ overload.

The lookup window must be at least the word size (k+1). The default window 7
works for common low-order models. Increasing the alphabet, order, or window
can substantially increase lookup-table and threshold-computation memory.

## Custom alphabets

Use `scan.Scanner` or `scan.scan` with an explicit alphabet. Each string lists
the characters mapped to a matrix row. The high-level `MotifScanner` uses the
DNA character mapping and is not the interface for arbitrary alphabets.

```python
from motifmatchpy import scan

alphabet = ["aA", "bB", "cC"]
matrix = [[2.0, -1.0], [-1.0, 2.0], [-1.0, -1.0]]
hits = scan.scan("ABCA", [matrix], [1 / 3] * 3, [4.0], 2, alphabet)[0]
assert [(hit.pos, hit.score) for hit in hits] == [(0, 4.0)]
```

Characters outside the alphabet split the sequence into scannable regions.
No biological complement mapping is inferred for custom alphabets.

## Explicit variants

Supply substitutions, insertions, or deletions through `Scanner.variant_matches`:

```python
from motifmatchpy import Variant, parsers, scan

matrix = parsers.pfm("docs/examples/toy.scores")
scanner = scan.Scanner(7)
scanner.set_motifs([matrix], [0.25] * 4, [6.0])
variants = [Variant(1, 2, "C")]
hits = scanner.variant_matches("ATG", variants, max_depth=1)[0]
assert [(hit.pos, hit.score, hit.variants) for hit in hits] == [(0, 6.0, (0,))]
```

`Variant(1, 1, "C")` inserts C; `Variant(1, 2, "")` deletes one base.
`max_depth` limits combined edits per match; zero means no explicit depth cap.
The method returns only matches enabled by edits. Call `scan()` separately for
ordinary hits. Read the [coordinate rules](formats.md#variant-coordinates)
before interpreting indel results.
