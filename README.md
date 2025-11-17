[![PyPI-Server](https://img.shields.io/pypi/v/biostrings.svg)](https://pypi.org/project/biostrings/)
![Unit tests](https://github.com/biocpy/biostrings/actions/workflows/run-tests.yml/badge.svg)

# biostrings

Efficient manipulation of genomic sequences in Python, inspired by the design of Bioconductor's [Biostrings](https://bioconductor.org/packages/Biostrings) package.

The core design relies on a **"pool and ranges"** memory model:

- **DNAStringSet** stores all sequences in a single contiguous block of memory (the pool).
- Individual sequences are defined by `start` and `width` coordinates (the ranges).
- Slicing a `DNAStringSet` returns a **view** (a new set of ranges pointing to the same pool), making subsetting operations virtually instantaneous and memory-free, regardless of the data size.

## Install

To get started, install the package from [PyPI](https://pypi.org/project/biostrings/)

```bash
pip install biostrings
```

## Quick Start

### Working with Single Sequences

The `DNAString` class represents a single DNA sequence. It enforces the IUPAC DNA alphabet and supports efficient byte-level operations.

```py
from biostrings import DNAString

# Create a DNA string
dna = DnaString("TTGAAAA-CTC-N")
print(dna)
# Output: TTGAAAA-CTC-N

# Basic operations
print(len(dna))            # 13
print(dna[0:3])            # DnaString(length=3, sequence='TTG')

# Reverse Complement
# Handles IUPAC ambiguity codes correctly (e.g., N -> N, M -> K)
rc = dna.reverse_complement()
print(rc)
# Output: N-GAG-TTTTCAA
```

### Working with Sets of Sequences

The `DNAStringSet` is the primary container for handling collections of sequences (e.g., reads from a FASTA file).

```py
from biostrings import DNAStringSet

# Efficiently create a set from a list of strings
seqs = [
    "ACGT",
    "GATTACA",
    "TTGAAAA-CTC-N",
    "ACGTACGT"
]
dss = DNAStringSet(seqs, names=["s1", "s2", "s3", "s4"])

print(dss)
# Output:
# <DNAStringSet of length 4>
#   [ 1]   4 ACGT                 s1
#   [ 2]   7 GATTACA              s2
#   [ 3]  13 TTGAAAA-CTC-N        s3
#   [ 4]   8 ACGTACGT             s4
```

<!-- biocsetup-notes -->

## Note

This project has been set up using [BiocSetup](https://github.com/biocpy/biocsetup)
and [PyScaffold](https://pyscaffold.org/).
