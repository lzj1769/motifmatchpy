# Python guide

Examples run from the repository root and use the [bundled files](quickstart.md).
The [API reference](reference/api.md) includes signatures and method descriptions.

## Load motifs and scan

```python
import motifmatchpy as mm

motif = mm.Motif.from_file("docs/examples/toy.scores", log_odds=False)
scanner = mm.MotifScanner([motif], threshold=6)
sequence = "ACGTTTACG"

for hit in scanner.scan(sequence, sequence_name="toy"):
    print(hit.name, hit.pos, hit.end, hit.strand, hit.score)
# toy.scores 0 3 + 6.0
# toy.scores 1 4 - 6.0
# toy.scores 6 9 + 6.0

assert scanner.count(sequence) == [3]
```

`Motif` is immutable and holds a name and scoring matrix. Its properties include
`order`, `width`, `length`, `min_score`, and `max_score`.
`motif.reverse_complement()` creates the opposite-strand motif.

`MotifScanner` searches both strands by default. Results are `Hit` objects sorted
by `(pos, motif.name, strand)`. `count()` returns counts in the order of the
input motifs, combining strands. `scan(..., max_hits=N)` keeps the first ordinary
hits per motif/strand rather than selecting the highest scores.

## Load several count matrices

```python
import motifmatchpy as mm

motifs = mm.read_motifs(["docs/examples/toy.pfm"], pseudocount=0.01)
scanner = mm.MotifScanner(motifs, p_value=0.02, bg="flat")
for hit in scanner.scan_fasta("docs/examples/toy.fa"):
    print(hit.sequence_name, hit.pos, hit.strand, hit.score)
```

`read_motifs` accepts multiple PFM/ADM paths, with optional `names` for explicit
motif labels. `scan_fasta()` reads one record at a time; it also accepts plain
text and gzip files. Each record's hits are collected before they are yielded.

## Configure thresholds

| Setting | Behavior |
| --- | --- |
| `p_value=1e-4` | Calculate a cutoff per motif and searched strand |
| `threshold=6.0` | Apply the same absolute cutoff to every scoring matrix |
| Neither cutoff | Only `scan_best_hits` is available |
| `bg="flat"` | Fixed uniform threshold background; Python default |
| `bg=[0.3, 0.2, 0.2, 0.3]` | Fixed threshold background |
| `bg="auto"` with `p_value` | Recalculate background and thresholds for each sequence |
| `both_strands=False` | Search only the supplied motif orientation |
| `window_size=7` | Lookup window; larger windows consume more memory |
| `threshold_precision=2000` | Score discretization for threshold calculation |

Do not set both `p_value` and `threshold`. Python defaults to a fixed uniform
background; the CLI's `scan -p` defaults to per-record estimation. To align them,
use `bg="auto"` in Python or `--batch` in the CLI.

`scanner.thresholds()` reports forward-motif thresholds. With `bg="auto"`, pass
a sequence, for example `scanner.thresholds("ACGTTTACG")`.
The conversion background is specified separately when loading a `Motif`.

## Approximate best hits

```python
import motifmatchpy as mm

motif = mm.Motif.from_file("docs/examples/toy.scores", log_odds=False)
scanner = mm.MotifScanner([motif])
hits = scanner.scan_best_hits("ACGTTTACG", target=1)
assert {(hit.pos, hit.strand) for hit in hits} == {(0, "+"), (1, "-"), (6, "+")}
```

Targets are approximate and applied independently to each motif and strand.
`iterations`, `mult`, and `limit_mult` tune the threshold search. This method
does not add IUPAC variant matches.

## IUPAC variants

```python
import motifmatchpy as mm

motif = mm.Motif.from_file("docs/examples/toy.scores", log_odds=False)
scanner = mm.MotifScanner([motif], threshold=6, both_strands=False)
hits = scanner.scan("AWG", include_variants=True)
assert hits == []  # W means A or T, so it cannot produce ACG.

hits = scanner.scan("AYG", include_variants=True)
assert len(hits) == 1  # Y means C or T; the C allele produces ACG.
assert hits[0].variant_sequence("AYG") == "aCg"
```

Ordinary hits have an empty `variants` tuple. Variant hits contain `Variant`
objects recording the substitutions used. `matched_sequence()` returns the
slice of the supplied sequence; `variant_sequence()` renders substitutions in
uppercase on a lowercase background. It does not reconstruct indels.

## Performance

Reuse a scanner when motifs and a fixed background stay unchanged. Setting up
lookup tables can dominate short scans. Auto backgrounds require rebuilding
thresholds and tables per record. Separate cold-start compilation time from
warm-run scanning time when benchmarking. Use `count()` when only counts are
needed, and avoid permissive cutoffs that generate very large hit lists.
