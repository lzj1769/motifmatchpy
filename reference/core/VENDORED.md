# Vendored MOODS core

Source: [MOODS 1.9.4.1](https://github.com/jhkorhonen/MOODS), `core/`, by Pasi
Rastas, Janne H. Korhonen and Petri Martinmäki. Dual-licensed under GPLv3 and
the Biopython license; see `COPYING.GPLv3` and `COPYING.BIOPYTHON` at the repo
root.

The SWIG interface files (`core/*.i`) are not vendored -- `../bindings.cpp`
replaces them with nanobind bindings.

`motifmatchpy` is a pure-Python reimplementation of these algorithms. This copy
exists only as a **test oracle**, so `tests/test_moods_parity.py` can show the
port is faithful on real data rather than only on hand-checked examples. It is
never a runtime dependency and is not published.

## Build settings

The extension is compiled with floating-point contraction disabled
(`-ffp-contract=off`, `/fp:precise` on MSVC). By default a compiler may fuse
`a + b * c` into a single multiply-add, which is more accurate but makes the
result depend on the compiler and the target CPU -- and puts it about one ULP
away from what Python computes. Turning contraction off is what lets the parity
tests assert bitwise equality rather than "close enough".

## Patches

Kept as small as possible so the copy stays easy to re-sync. Each is marked in
the source with a `motifmatchpy patch:` comment, and each has a corresponding
regression test.

### 1. `motif.h` -- missing virtual destructor

`Motif` is an abstract base whose instances are owned through
`std::unique_ptr<Motif>` in `Scanner::motifs`. With no virtual destructor,
destroying one is undefined behaviour and leaks the derived class's members
(score matrix, lookahead tables) on every `set_motifs()` call.

### 2. `scanner.h` -- `Scanner` declared copyable but is not

`Scanner` holds a `vector<unique_ptr<Motif>>`, so copying it never compiles --
yet `std::is_copy_constructible<Scanner>` reports true, because the vector's
copy constructor is only ill-formed once instantiated. Generic code that trusts
the trait (nanobind's, among others) then fails to build. The copy operations
are now deleted explicitly and the move operations defaulted.

### 3. `moods_parsers.cpp` -- out-of-bounds read on an unparseable file

`pfm()` and `pfm_to_log_odds()` evaluate `mat[0].size()` before checking whether
the table has any rows, so a missing, empty or unreadable matrix file indexes
past the end of an empty vector instead of reporting a parse failure. The
emptiness check now runs first.

### 4. `moods_parsers.cpp` -- `adm_to_log_odds` ignored its log base

The function takes `log_base` and then never uses it, always applying the
natural logarithm. `--log-base` therefore rescaled `.pfm` matrices but not
`.adm` ones, putting the two on different scales -- and their scores and
thresholds out of any common frame -- within a single run.

### 5. `moods_tools.cpp` -- wrong symbol read in the high-order p-value DP

Setting up the score distribution for a high-order matrix, the background
weight of a q-gram was computed as

```cpp
prob *= bg[(CODE >> (q - i - 1)) & A_MASK];      // upstream
prob *= bg[(CODE >> (SHIFT * (q - i - 1))) & A_MASK];  // fixed
```

The shift is in bits, not symbols, so upstream reads the wrong symbol out of the
code (its own `// TODO: check correctness` sits three lines above). The weights
only happen to sum to 1 when the background is uniform; with
`bg = [0.7, 0.1, 0.1, 0.1]` they sum to 1.36, and every threshold derived from a
non-uniform background is wrong.

### 6. `scanner.cpp` -- matches dropped at the end of a scannable region

In `process_matches`, the loop handling positions near the end of a region
guarded a match with `size() < end - i`, while the equivalent guard for very
short sequences a few lines above uses `<=`. A motif occupies `[i, i + size)`
and so fits exactly when `i + size <= end`, making `<` one position too strict.

For any motif **shorter than the scanning window** (default 7), a match ending
exactly at the end of a scannable region was silently discarded -- the end of
the sequence, or the base before an `N`. Real genomic sequence is full of `N`
runs, so every one of them was an opportunity to lose a match.

```
sequence           window  upstream  fixed
CCCCCCCCCCAAAA          7  []        [10]
CCCCCCAAAANCCCCCC       7  []        [6]
```

### 7. `scanner.cpp` -- matches lost in regions shorter than the window

The branch for a region shorter than the window tested the lookup table with the
code from the *previous* position and only then advanced the window:

```cpp
if (match_handler.has_hits(code))       // previous position's window
{
    code = (code << SHIFT) & MASK;      // ... then move to this one
    for (const scanner_output& y : match_handler.hits(code))
```

So each position's lookup was gated on its predecessor's, and the window stood
still whenever the predecessor had no candidates. Any match not starting at the
first position of such a region was lost: with a 3 bp motif and window 7,
`CAAA` found nothing. Advancing first and then testing fixes it; an exhaustive
sweep of all 349,524 sequences up to length 9 now agrees with a brute-force scan
at every window size.

### 8. `motif_h.cpp` -- crashes when the window does not fit the motif

Two separate faults for a high-order motif shorter than the scanning window,
which is what the default window of 7 gives you for any short adjacent
dinucleotide model:

* `max_scores_f`/`max_scores_b` hold the window size in a `double` and pass it
  to a vector constructor. When `end < start` that is a negative value converted
  to `size_t`, so `l > m` raised `std::length_error` instead of producing an
  empty table.
* A window smaller than the motif's q-gram length left `window_position()`
  calling `.back()` on an empty table -- undefined behaviour, in practice a
  segfault. It is a meaningless configuration, so the constructor now rejects it
  with `std::invalid_argument`.

### 9. `moods_scan.cpp` -- high-order `naive_scan_dna` stopped one position early

Its loop ran while `i + cols + q - 1 < end`, where the 0-order overload uses
`i + m <= end`, so a motif ending exactly at the end of the sequence was missed.
This affects only the brute-force reference helper, but it made that helper
disagree with the scanner it exists to check.

### 10. `motif_0.cpp` -- lookahead order left to an unstable sort

`compute_lookahead_order` ranks the positions outside the window by how well
each discriminates, using `std::sort` with

```cpp
bool operator() (int i, int j) { return (*ed)[i] > (*ed)[j]; }
```

Positions that discriminate equally well compare equal, which the JASPAR
matrices produce often -- 36 of the matrix, background and window-size
combinations the tests cover have at least one tie inside the lookahead set.
`std::sort` is not stable, so tied positions come out in whatever order the
implementation happens to produce, and `check_hit` adds the scores up in exactly
that order. The last bits of every score therefore depended on which standard
library built the binary: libc++ and libstdc++ disagree, and the parity tests
passed on macOS while failing on Linux.

`std::stable_sort` makes the order the same everywhere, and the same as the
(always stable) Python sort. Scores shift by at most a ULP -- 15 to 159 of a few
hundred matches, in the tests that exposed it -- and no match position moves.

## Known upstream behaviour left as-is

* **256-entry symbol tables are indexed with a plain `char`.** `snp_variants`,
  `bg_from_sequence_dna`, `Scanner::process_matches` and both `check_hit`
  implementations index tables of size 256 with `seq[i]`, which is signed on the
  usual platforms, so any byte above 127 reads outside the table. Rather than
  patch every site, both `motifmatchpy` and this wrapper reject non-ASCII
  sequences at the Python boundary.

* **`MotifH::window_position` assigns where it means to accumulate.** Its
  running expectation is built with `current_exp = es[i]` inside the loop over
  the initial window, so it ends up holding one column's mean rather than their
  sum. The choice of window only steers the search, never which positions match,
  so this is reproduced rather than fixed: matching it keeps the order in which
  scores are added up the same, and with it the last bits of every score.
