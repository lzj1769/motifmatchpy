# Sequence I/O

`read_fasta` and `read_sequences` are also exported from `motifmatchpy`.
They yield `(name, sequence)` pairs. Compression is detected from file contents.

::: motifmatchpy.io
    options:
      members: [read_fasta, read_sequences, open_text]
