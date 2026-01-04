from __future__ import annotations

import re
from copy import deepcopy
from typing import Any, Dict, Optional, Union
from warnings import warn

import biocutils as ut

from .utils import _sanitize_metadata

# From R's DNA_ALPHABET
# We'll use the standard IUPAC DNA codes + gap
DNA_IUPAC_LETTERS = "ACGTRYSWKMBDHVN-"
DNA_IUPAC_BYTES = b"ACGTRYSWKMBDHVN-"

# Pre-compiled regex for validation
_DNA_VALIDATOR = re.compile(f"^[{''.join(DNA_IUPAC_LETTERS)}]*$", re.IGNORECASE)

# Translation tables for reverse_complement
_DNA_COMPLEMENT_TABLE = str.maketrans("ACGTRYSWKMBDHVN-", "TGCAYRSWMKVHDBN-")
_DNA_COMPLEMENT_TABLE_BYTES = bytes.maketrans(b"ACGTRYSWKMBDHVN-acgtryswkmbdhvn-", b"TGCAYRSWMKVHDBN-TGCAYRSWMKVHDBN-")

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


class DNAString(ut.BiocObject):
    """A string container for a DNA sequence, similar to Bioconductor's DNAString.

    This class stores the sequence internally as bytes, enforcing the
    DNA alphabet.
    """

    def __init__(
        self,
        sequence: Union[str, bytes],
        metadata: Optional[Union[Dict[str, Any], ut.NamedList]] = None,
        _validate: bool = True,
    ):
        """Create a DNAString.

        Args:
            sequence:
                A string or bytes object representing a DNA sequence.

            metadata:
                Additional metadata. If None, defaults to an empty dictionary.

            _validate:
                Whether to validate the arguments, internal use only.
        """
        super().__init__(metadata=metadata, _validate=_validate)

        if isinstance(sequence, str):
            self._data = sequence.upper().encode("ascii")
        elif isinstance(sequence, bytes):
            self._data = sequence
        else:
            raise TypeError(f"Cannot initialize DNAString with type {type(sequence)}")

        self._metadata = _sanitize_metadata(metadata)

        if _validate:
            if not _DNA_VALIDATOR.match(self._data.decode("ascii")):
                raise ValueError("Input string contains non-DNA characters.")

    #################
    #### Copying ####
    #################

    def __copy__(self) -> DNAString:
        """Shallow copy of the object.

        Returns:
            Same type as the caller, a shallow copy of this object.
        """
        return type(self)(
            sequence=str(self),
            metadata=self._metadata,
            _validate=False,
        )

    def __deepcopy__(self, memo) -> DNAString:
        """Deep copy of the object.

        Args:
            memo: Passed to internal :py:meth:`~deepcopy` calls.

        Returns:
            Same type as the caller, a deep copy of this object.
        """
        return type(self)(
            sequence=deepcopy(str(self), memo),
            metadata=deepcopy(self._metadata, memo),
            _validate=False,
        )

    ########################
    #### Getter/setters ####
    ########################

    def get_sequence(self) -> str:
        """Get the sequence.

        Returns:
            The sequence string.
        """
        return self._data.decode("ascii")

    ##################
    #### printing ####
    ##################

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
        return f"DNAString(length={length}, sequence='{snippet}')"

    #####################
    #### comparators ####
    #####################

    def __len__(self) -> int:
        """Return the length of the sequence."""
        return len(self._data)

    def __eq__(self, other) -> bool:
        """Check for equality with another DNAString or str."""
        if isinstance(other, DNAString):
            return self._data == other._data
        if isinstance(other, str):
            return str(self) == other.upper()
        return False

    #########################
    #### Getitem/setitem ####
    #########################

    def __getitem__(self, key: Union[int, slice]) -> DNAString:
        """Extract a subsequence (slicing).

        Args:
            key:
                An integer or slice.

        Returns:
            A new DNAString object representing the subsequence.
        """
        if isinstance(key, int):
            if key < 0:
                key += len(self)
            return DNAString(self._data[key : key + 1])
        elif isinstance(key, slice):
            return DNAString(self._data[key])
        else:
            raise TypeError(f"Index must be int or slice, not {type(key)}.")

    #################
    #### methods ####
    #################

    def reverse_complement(self) -> DNAString:
        """Compute the reverse complement of the sequence.

        Returns:
            A new DNAString with the reverse complement.
        """
        complemented = self._data.translate(_DNA_COMPLEMENT_TABLE_BYTES)
        return DNAString(complemented[::-1])

    def to_bytes(self) -> bytes:
        """Get the underlying byte representation."""
        return self._data
