# Scan ATAC-seq peaks against a reference genome

Provide a peak BED file, the matching genome assembly's FASTA, and motif
files in JASPAR, MEME, TRANSFAC, PFM or ADM format. `scan` extracts each interval, scans both strands, and writes one
BED6 row per motif hit with **genomic coordinates**.

## Typical command

```text
motifmatchpy scan --regions atac_peaks.bed --genome hg38.fa \
  -m JASPAR2026_CORE_vertebrates.jaspar -p 1e-4 -o motif_hits.bed
```

No intermediate FASTA extraction step, organism database installation, or genome
index is needed. The FASTA is read once, retaining one requested chromosome at
a time. Memory therefore scales with the largest requested chromosome, not just
the peak length. Plain and gzip-compressed inputs are supported. Unrequested
chromosomes are read past without retaining their sequence; an existing `.fai`
is not used.

Use the same assembly for peaks and reference, and exactly matching chromosome
identifiers (`chr1` and `1` are different). Missing chromosomes, invalid BED
coordinates, and intervals beyond a chromosome's end produce errors rather
than silently clipping or skipping peaks. On an error, a partially written
output file may remain; check the command's exit status before using it.

## Run the bundled example

From the repository root:

```sh
motifmatchpy scan --regions docs/examples/peaks.bed \
  --genome docs/examples/genome.fa -m docs/examples/motifs.jaspar -p 0.02
```

```text
chr1	2	5	MA0000.1_TEST	4.15663	+
chr1	3	6	MA0000.1_TEST	4.15663	-
chr1	8	11	MA0000.1_TEST	4.15663	+
```

The fixture contains two motifs; only `TEST` matches these peaks. The input
peak `[2, 6)` contains a forward `ACG` and reverse `CGT`. The last hit comes
from `[8, 11)`. Coordinates are zero-based and half-open on both strands.

## Other motif file formats

The same command works with MEME or TRANSFAC input; no format flag is needed:

```sh
motifmatchpy scan --regions docs/examples/peaks.bed \
  --genome docs/examples/genome.fa -m docs/examples/motifs.meme -p 0.02
motifmatchpy scan --regions docs/examples/peaks.bed \
  --genome docs/examples/genome.fa -m docs/examples/motifs.transfac -p 0.02
```

You can also pass several files after `-m`, even when they use different formats.
All records are loaded. See [file formats](formats.md) for MEME probability
conversion and naming rules.

## Parameters and defaults

| Option | Purpose |
| --- | --- |
| `--regions BED` | BED3 or additional-column BED, including narrowPeak |
| `--genome FASTA` | Reference sequences; FASTA identifier is the first header token |
| `-m`, `--matrices FILE ...` | Auto-detected JASPAR, MEME, TRANSFAC, numeric PFM or ADM |
| `-S FILE ...` | Existing score matrices, also usable with peaks |
| `-p P` | P-value cutoff; default `1e-4` when no cutoff is supplied in peak mode |
| `-t T` / `-B N` | Alternatives: absolute threshold / approximate best hits per peak |
| `--bg A C G T` | Fixed threshold background; default uniform |
| `--lo-bg A C G T` | Count-conversion background; default uniform |
| `--ps P` | Conversion pseudocount; default `0.01` |
| `-R` | Forward-strand search only |
| `-o FILE` | Output file; default stdout |
| `-f bed\|csv\|tsv` | Default BED in peak mode |

`--regions` and `--genome` must be given together and cannot be combined with `-s`.
Peak mode always uses a fixed background, so `--batch` is implicit. The same
prepared scanner is reused across peaks. Original `-s` mode retains its
per-record background default and requires an explicit cutoff.

## Output and overlapping peaks

BED6 fields are `chrom`, `start`, `end`, `motif`, `score`, `strand`. JASPAR names
combine header ID and TF name, such as `MA0139.1_CTCF`. Scores are raw log-odds
values, not normalized 0–1000 BED scores. The default logarithm is natural log;
use `--log-base 2` for bits. CSV/TSV also uses genomic positions in peak mode,
while its `match` column shows the actual matched reference bases.

Only hits entirely contained within an individual peak are reported. Peaks are
not joined, extended, or centered on summits. Additional BED columns, including
peak name and strand, do not change scanning. Each peak is scanned independently,
so identical or overlapping peaks can produce duplicate hits. Both strands of
a palindromic motif are retained.

Records follow FASTA chromosome order, then peak start/end order, with hits
sorted within each peak. Overlapping peaks mean the resulting BED is not
necessarily globally sorted. To obtain sorted, unique rows on Unix:

```sh
sort -k1,1 -k2,2n -k3,3n motif_hits.bed | uniq > motif_hits.unique.bed
```

For peak-associated analyses, retain duplicates or preserve the original peaks
and intersect hits back to them, rather than discarding associations blindly.

## Relation to RGT matching

This workflow follows the region-to-motif-instance BED pattern of
[`rgt-motifanalysis matching`](https://reg-gen.readthedocs.io/en/latest/motif_analysis/tool_usage.html).
Here the genome and motif collection are supplied directly, rather than through
an organism installation. It does not implement RGT enrichment, random regions,
BigBed conversion, promoter selection, or all of RGT's filtering options.
Defaults such as pseudocounts differ; do not assume numeric equivalence without
aligning matrix conversion, backgrounds, and threshold settings.
