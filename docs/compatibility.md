# Compatibility and troubleshooting

## Relationship to MOODS

motifmatchpy implements the main MOODS 1.9.4.1 C++ algorithm families, including
multi-matrix DNA scanning, higher-order scoring, p-value thresholds, approximate
best hits, and variant matching. It adds a high-level Python API and the
`scan`, `threshold`, and `info` subcommands.

The bundled `reference/` is a **patched**, development-only C++ test oracle. It
is not a runtime dependency or an unmodified upstream binary. Corrections
include boundary matches, high-order threshold handling, and explicit ADM log
bases. The exact changes are recorded in
[VENDORED.md](https://github.com/lzj1769/motifmatchpy/blob/main/reference/core/VENDORED.md).
Do not assume universal bit-for-bit equality with unmodified upstream binaries.

The low-level modules have familiar MOODS names, but this is not a drop-in
replacement for every SWIG type, overloaded signature, or `moods-dna.py` option.
Use the documented Python signatures and CLI options. In particular, CLI
best-hits does not chain into SNP scanning; use a fixed cutoff for variants.

## No hits

- Check `info`: is the threshold above the motif's maximum score?
- Check whether you passed counts as `-S` or scores as `-m` accidentally.
- Try the [tiny example](quickstart.md) to verify installation.
- Check sequence length, unknown bases, and `-R` strand filtering.
- Remember that very small p-values can exclude every position.

## Different results between CLI and Python

Align the matrix conversion background, threshold background, pseudocount,
logarithm base, p-value precision, strands, and IUPAC behavior. CLI p-value scans
estimate background per record by default; Python uses a fixed uniform
background. CLI cutoff scans include IUPAC alternatives by default; Python's
`scan()` requires `include_variants=True`.

CSV/TSV scores have six significant digits, whereas Python retains full
floating-point values. Positions are zero-based, and reverse-strand match text
stays in the input orientation.

## Slow first scan or high memory usage

The first scan can compile Numba kernels. Reuse the scanner and compare warm
runs separately from compilation. Larger lookup windows, high-order matrices,
more permissive cutoffs, and larger threshold precision can consume more
memory. Start with the default window and precision, and use a fixed background
when appropriate.

## Matrix parsing errors

PFM input must have four numeric rows for DNA. Row labels and brackets are not
supported. ADM input needs 16 conditional rows plus four nonempty initial-term
rows. Higher-order numeric score matrices should be constructed through the
Python API rather than treated as ordinary PFM files.

## Reporting a problem

Open an [issue](https://github.com/lzj1769/motifmatchpy/issues) with the package
version, Python/platform details, exact command or small script, minimal matrix
and sequence inputs, and expected versus observed results. For upstream
comparisons, include the MOODS version, compiler settings, and whether the
reference contains patches.
