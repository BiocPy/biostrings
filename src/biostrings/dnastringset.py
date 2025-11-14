import numpy as np
from biocpy.iranges import IRanges

from .DnaString import _DNA_VALIDATOR, DnaString


class DnaStringset:
    """A collection of DNA sequences, similar to Bioconductor's DNAStringSet.

    This class follows the "pool and ranges" model for high memory
    efficiency. All sequences are stored in a single concatenated 'bytes' object (the pool).

    An 'IRanges' object tracks the start and width of each sequence in the pool.
    """

    def __init__(
        self,
        sequences: list[str],
        names: list[str] = None,
        _pool: bytes = None,
        _ranges: IRanges = None,
    ):
        """Create a DnaStringset.

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
        """
        if _pool is not None and _ranges is not None:
            # Internal constructor: just assign the shared pool and new ranges
            self._pool = _pool
            self._ranges = _ranges
        elif sequences is not None:
            # Public constructor: build the pool and ranges
            pool_parts = []
            widths = np.zeros(len(sequences), dtype=np.int32)

            for i, seq_str in enumerate(sequences):
                if not _DNA_VALIDATOR.match(seq_str):
                    raise ValueError(f"Sequence at index {i} contains non-DNA characters.")

                seq_bytes = seq_str.upper().encode("ascii")
                pool_parts.append(seq_bytes)
                widths[i] = len(seq_bytes)

            self._pool = b"".join(pool_parts)
            starts = np.concatenate((np.array([0], dtype=np.int32), np.cumsum(widths[:-1])))

            self._ranges = IRanges(starts=starts, widths=widths, names=names)
        else:
            # Empty set
            self._pool = b""
            self._ranges = IRanges([], [], names=[])

    def __len__(self) -> int:
        """Return the number of sequences in the set."""
        return len(self._ranges)

    def width(self) -> np.ndarray:
        """Return an array of lengths for all sequences."""
        return self._ranges.width

    @property
    def names(self) -> list[str]:
        """Return the names of the sequences."""
        return self._ranges.names

    @names.setter
    def names(self, new_names: list[str]):
        """Set the names of the sequences."""
        self._ranges.names = new_names

    def __getitem__(self, key: Union[int, slice, list[int], np.ndarray]) -> Union[DnaString, "DnaStringset"]:
        """
        Extract one or more sequences.

        - If key is int: Returns a DnaString object (a copy).
        - If key is slice or list: Returns a new DnaStringset (a view).
        """
        if isinstance(key, int):
            if key < 0:
                key += len(self)

            r = self._ranges.ranges[key]
            start = r["start"]
            end = start + r["width"]

            # Return a DnaString, which is a *copy* of the bytes
            return DnaString(self._pool[start:end])

        elif isinstance(key, (slice, list, np.ndarray)):
            # New IRanges object (subsetted)
            new_ranges = self._ranges[key]

            # Return a new DnaStringset *view*
            # It shares the *same* self._pool, but has new ranges.
            # This is the memory-efficient design.
            return DnaStringset(sequences=None, _pool=self._pool, _ranges=new_ranges)
        else:
            raise TypeError(f"Index must be int, slice, or list, not {type(key)}")

    def to_list(self) -> list[str]:
        """
        Convert the set to a list of Python strings.

        [See R: new_CHARACTER_from_XStringSet in XStringSet_class.c]
        This is a candidate for C++ optimization.
        """
        output = []
        for r in self._ranges.ranges:
            start = r["start"]
            end = start + r["width"]
            output.append(self._pool[start:end].decode("ascii"))
        return output

    def unlist(self) -> DnaString:
        """
        Concatenate all sequences in the set into one DnaString.

        [See R: XStringSet_unlist in XStringSet_class.c]
        """
        # This is fast if the pool is already concatenated, but
        # if the ranges are subsetted, we need to build it.
        if len(self) == 0:
            return DnaString("")

        # Check if the ranges cover the pool contiguously
        first_start = self._ranges.ranges[0]["start"]
        last_r = self._ranges.ranges[-1]
        last_end = last_r["start"] + last_r["width"]

        if last_end - first_start == np.sum(self.width()):
            # Simple case: just a view on the pool
            return DnaString(self._pool[first_start:last_end])
        else:
            # Complex case: ranges are subsetted, must join
            return DnaString(b"".join(self.to_list()))

    def __repr__(self) -> str:
        """Return a compact representation."""
        n = len(self)
        cls_name = self.__class__.__name__

        if n == 0:
            return f"<{cls_name} of length 0>"

        header = f"<{cls_name} of length {n}>"

        # Show logic from R
        max_show = 10
        half_show = 5

        lines = []
        widths = self.width()
        names = self.names if self.names else [""] * n
        max_width_str = len(str(np.max(widths)))

        def format_line(i):
            w = widths[i]
            seq = self[i]  # This gets a DnaString

            if w > 18:
                snippet = str(seq[:7]) + "..." + str(seq[-8:])
            else:
                snippet = str(seq)

            name_str = names[i]
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
