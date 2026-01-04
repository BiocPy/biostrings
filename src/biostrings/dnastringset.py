from __future__ import annotations

from copy import deepcopy
from typing import Any, Dict, List, Optional, Union
from warnings import warn

import biocutils as ut
import numpy as np
from iranges import IRanges

from .dnastring import DNAString
from .utils import _sanitize_metadata

try:
    from . import lib_biostrings

    _CPP_OPS_ENABLED = True
except ImportError:
    _CPP_OPS_ENABLED = False
    from .dnastring import _DNA_VALIDATOR


__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


class DNAStringSet(ut.BiocObject):
    """A collection of DNA sequences, similar to Bioconductor's DNAStringSet.

    This class follows the "pool and ranges" model for high memory
    efficiency. All sequences are stored in a single concatenated 'bytes' object (the pool).

    An 'IRanges' object tracks the start and width of each sequence in the pool.
    """

    def __init__(
        self,
        sequences: Optional[List[str]] = None,
        names: Optional[Union[List[str], ut.Names]] = None,
        _pool: Optional[bytes] = None,
        _ranges: Optional[IRanges] = None,
        metadata: Optional[Union[Dict[str, Any], ut.NamedList]] = None,
        _validate: bool = True,
    ):
        """Create a DNAStringSet.

        Args:
            sequences:
                A list of Python strings to initialize the set.

            names:
                An optional list of names for the sequences.

            _pool (internal):
                Used by methods like __getitem__ to create
                new sets without copying data.

            _ranges (internal):
                Used by methods like __getitem__.

            metadata:
                Additional metadata. If None, defaults to an empty dictionary.

            validate:
                Whether to validate the arguments, internal use only.
        """
        super().__init__(metadata=metadata, _validate=_validate)

        if _pool is not None and _ranges is not None:
            self._pool = _pool
            self._ranges = _ranges

        elif sequences is not None and len(sequences) > 0:
            if _validate:
                for i, seq_str in enumerate(sequences):
                    if not _DNA_VALIDATOR.match(seq_str):
                        raise ValueError(f"Sequence at index {i} contains non-DNA characters.")

            if _CPP_OPS_ENABLED:
                try:
                    self._pool, starts, widths = lib_biostrings.create_DNAStringSet_pool(sequences)
                    self._ranges = IRanges(starts=starts, widths=widths, names=names)
                except Exception as e:
                    raise e
            else:
                pool_parts = []
                widths = np.zeros(len(sequences), dtype=np.int32)

                for i, seq_str in enumerate(sequences):
                    seq_bytes = seq_str.upper().encode("ascii")
                    pool_parts.append(seq_bytes)
                    widths[i] = len(seq_bytes)

                self._pool = b"".join(pool_parts)
                starts = np.concatenate((np.array([0], dtype=np.int32), np.cumsum(widths[:-1])))

                self._ranges = IRanges(start=starts, width=widths, names=names)

        else:
            # Empty set
            self._pool = b""
            self._ranges = IRanges(np.array([], dtype=np.int32), np.array([], dtype=np.int32), names=[])

        self._metadata = _sanitize_metadata(metadata)

    #################
    #### Copying ####
    #################

    def __copy__(self) -> DNAStringSet:
        """Shallow copy of the object.

        Returns:
            Same type as the caller, a shallow copy of this object.
        """
        return type(self)(
            _pool=self._pool,
            _ranges=self._ranges,
            metadata=self._metadata,
            _validate=False,
        )

    def __deepcopy__(self, memo) -> DNAStringSet:
        """Deep copy of the object.

        Args:
            memo: Passed to internal :py:meth:`~deepcopy` calls.

        Returns:
            Same type as the caller, a deep copy of this object.
        """
        return type(self)(
            _pool=deepcopy(self._pool, memo),
            _ranges=deepcopy(self._ranges, memo),
            metadata=deepcopy(self._metadata, memo),
            _validate=False,
        )

    ########################
    #### Getter/setters ####
    ########################

    def get_names(self) -> Optional[ut.Names]:
        """Get range names.

        Returns:
            List containing the names for all ranges, or None if no names are
            present.
        """
        return self._ranges.get_names()

    def set_names(self, names: Optional[List[str]], in_place: bool = False) -> DNAStringSet:
        """
        Args:
            names:
                Sequence of names or None, see the constructor for details.

            in_place:
                Whether to modify the object in place.

        Returns:
            If ``in_place = False``, a new ``DNAStringSet`` is returned with the
            modified names. Otherwise, the current object is directly modified
            and a reference to it is returned.
        """
        output = self._define_output(in_place)
        output._ranges.set_names(names, in_place=in_place)
        return output

    @property
    def names(self) -> Optional[ut.Names]:
        """Return the names of the sequences."""
        return self._ranges.get_names()

    @names.setter
    def names(self, new_names: List[str]):
        """Set the names of the sequences."""
        warn(
            "Setting property 'names' is an in-place operation, use 'set_names' instead",
            UserWarning,
        )
        self.set_names(names=new_names, in_place=True)

    ##################
    #### printing ####
    ##################

    def __repr__(self) -> str:
        """Return a compact representation."""
        n = len(self)
        cls_name = self.__class__.__name__

        if n == 0:
            return f"<{cls_name} of length 0>"

        header = f"<{cls_name} of length {n}>"

        max_show = 10
        half_show = 5

        lines = []
        widths = self.width()
        names = self.names if self.names else [""] * n
        max_width_str = len(str(np.max(widths))) if len(widths) > 0 else 0

        def format_line(i):
            w = widths[i]
            # TODO: Avoid creating DnaString object just for repr
            # Access bytes directly
            start = self._ranges.start[i]
            end = start + w
            seq_bytes = self._pool[start:end]

            if w > 18:
                snippet = seq_bytes[:7].decode("ascii") + "..." + seq_bytes[-8:].decode("ascii")
            else:
                snippet = seq_bytes.decode("ascii")

            name_str = names[i] if names[i] else ""
            if len(name_str) > 10:
                name_str = name_str[:7] + "..."

            return f"  [{i + 1:2d}] {w:>{max_width_str}d} {snippet:<20} {name_str}"

        if n <= max_show:
            for i in range(n):
                lines.append(format_line(i))
        else:
            for i in range(half_show):
                lines.append(format_line(i))
            lines.append(f"  ... {n - (2 * half_show)} more sequences ...")
            for i in range(n - half_show, n):
                lines.append(format_line(i))

        return header + "\n" + "\n".join(lines)

    ##################
    #### the rest ####
    ##################

    def __len__(self) -> int:
        """Return the number of sequences in the set."""
        return len(self._ranges)

    def get_width(self) -> np.ndarray:
        """Return an array of lengths for all sequences."""
        return self._ranges.get_width()

    def width(self) -> np.ndarray:
        """Alias to :py:meth:`~.get_width`."""
        return self.get_width()

    def __getitem__(self, key: Union[int, slice, List[int], np.ndarray]) -> Union[DNAString, DNAStringSet]:
        """Extract one or more sequences.

        Args:
            key:
                - If key is int: Returns a DNAString object (a copy).
                - If key is slice or list: Returns a new DNAStringSet (a view).

        Returns:
            A DNAString or DNAStringSet object representing the slice.
        """
        if isinstance(key, int):
            if key < 0:
                key += len(self)

            r = self._ranges[key]
            start = r._start[0]
            end = start + r._width[0]

            return DNAString(self._pool[start:end])

        elif isinstance(key, (slice, list, np.ndarray)):
            new_ranges = self._ranges[key]
            return DNAStringSet(sequences=None, _pool=self._pool, _ranges=new_ranges)
        else:
            raise TypeError(f"Index must be int, slice, or list, not {type(key)}")

    def to_list(self) -> List[str]:
        """Convert the set to a list of Python strings."""
        output = []
        for i in range(len(self._ranges)):
            start = self._ranges._start[i]
            width = self._ranges._width[i]
            end = start + width
            output.append(self._pool[start:end].decode("ascii"))
        return output

    def unlist(self) -> DNAString:
        """Concatenate all sequences in the set into one DNAString."""
        if len(self) == 0:
            return DNAString("")

        first_start = self._ranges._start[0]
        last_start = self._ranges._start[-1]
        last_width = self._ranges._width[-1]
        last_end = last_start + last_width

        total_width = np.sum(self.width())

        if last_end - first_start == total_width:
            return DNAString(self._pool[first_start:last_end])
        else:
            return DNAString(b"".join([x.encode("ascii") for x in self.to_list()]))
