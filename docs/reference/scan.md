# Scanning

The low-level API returns one group per input matrix. It does not automatically
add reverse complements. [Examples](../advanced.md) cover DNA, custom alphabets,
and explicit variants. `Scanner` and the DNA helpers are also exported at the
package root.

`scan_max_hits` keeps the first hits up to the limit; it does not rank by score.
Variant coordinates refer to the edited sequence.

::: motifmatchpy.scan
    options:
      members: [Scanner, scan_dna, scan, scan_best_hits_dna, naive_scan_dna]
