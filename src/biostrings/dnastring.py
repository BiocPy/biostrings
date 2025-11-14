import re
from typing import Union

# From R's DNA_ALPHABET
# We'll use the standard IUPAC DNA codes + gap
DNA_IUPAC_LETTERS = "ACGTRYSWKMBDHVN-"
DNA_IUPAC_BYTES = b"ACGTRYSWKMBDHVN-"

# Pre-compiled regex for validation
_DNA_VALIDATOR = re.compile(f"^[{''.join(DNA_IUPAC_LETTERS)}]*$", re.IGNORECASE)

# Translation tables for reverse_complement
_DNA_COMPLEMENT_TABLE = str.maketrans("ACGTRYSWKMBDHVN-", "TGCAYRSWMKVHDBN-")
_DNA_COMPLEMENT_TABLE_BYTES = bytes.maketrans(b"ACGTRYSWKMBDHVN-acgtryswkmbdhvn-", b"TGCAYRSWMKVHDBN-TGCAYRSWMKVHDBN-")


class DnaString:
    """A string container for a DNA sequence, similar to Bioconductor's DNAString.

    This class stores the sequence internally as bytes, enforcing the
    DNA alphabet.
    """

    def __init__(self, sequence: Union[str, bytes]):
        """Create a DnaString.

        Args:
            sequence:
                A string or bytes object representing a DNA sequence.
        """
        if isinstance(sequence, str):
            # Validate string input
            if not _DNA_VALIDATOR.match(sequence):
                raise ValueError("Input string contains non-DNA characters.")
            self._data = sequence.upper().encode("ascii")
        elif isinstance(sequence, bytes):
            # Assume bytes are already valid (e.g., from internal slicing)
            self._data = sequence
        else:
            raise TypeError(f"Cannot initialize DnaString with type {type(sequence)}")

    def __len__(self) -> int:
        """Return the length of the sequence."""
        return len(self._data)

    def __str__(self) -> str:
        """Return the sequence as a Python string."""
        return self._data.decode("ascii")

    def __repr__(self) -> str:
        """Return a string representation."""
        length = len(self)
        if length > 20:
            snippet = str(self[:10]) + "..." + str(self[-10:])
        else:
            snippet = str(self)
        return f"DnaString(length={length}, sequence='{snippet}')"

    def __eq__(self, other) -> bool:
        """Check for equality with another DnaString or str."""
        if isinstance(other, DnaString):
            return self._data == other._data
        if isinstance(other, str):
            return str(self) == other.upper()
        return False

    def __getitem__(self, key: Union[int, slice]) -> "DnaString":
        """Extract a subsequence (slicing).

        Args:
            key:
                An integer or slice.

        Returns:
            A new DnaString object representing the subsequence.
        """
        if isinstance(key, int):
            if key < 0:
                key += len(self)
            return DnaString(self._data[key : key + 1])
        elif isinstance(key, slice):
            return DnaString(self._data[key])
        else:
            raise TypeError(f"Index must be int or slice, not {type(key)}")

    def reverse_complement(self) -> "DnaString":
        """Compute the reverse complement of the sequence.

        Returns:
            A new DnaString with the reverse complement.
        """
        complemented = self._data.translate(_DNA_COMPLEMENT_TABLE_BYTES)
        return DnaString(complemented[::-1])

    def to_bytes(self) -> bytes:
        """Get the underlying byte representation."""
        return self._data
