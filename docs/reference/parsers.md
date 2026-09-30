# Matrix parsers

These functions read numeric matrices, ADM terms, and JASPAR/MEME/TRANSFAC
collections. `jaspar()`, `meme()`, and `transfac()` return `(name, count_matrix)`
pairs. `read_motifs()` auto-detects formats and converts counts to scores. See [file formats](../formats.md).

::: motifmatchpy.parsers
    options:
      members: [detect_format, jaspar, meme, transfac, pfm, pfm_to_log_odds, adm_1o_terms, adm_0o_terms, adm_to_log_odds, read_table, ParseError]
