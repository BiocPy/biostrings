import numpy as np
import pytest

from biostrings import DNAString, DNAStringSet


@pytest.fixture
def sample_seqs():
    """Provides a sample set of sequences and names."""
    return {
        "seqs": [
            "ACGT",
            "GATTACA",
            "",
            "TTGAAAA-CTC-N",  #
            "ACGTACGT",
        ],
        "names": ["seq1", "seq2", "empty", "iupac", "seq5"],
    }


@pytest.fixture
def sample_set(sample_seqs):
    """Creates a DNAStringSet from the sample_seqs fixture."""
    return DNAStringSet(sample_seqs["seqs"], names=sample_seqs["names"])


def test_construction(sample_set, sample_seqs):
    assert len(sample_set) == 5
    assert list(sample_set.names) == sample_seqs["names"]
    assert np.array_equal(sample_set.width(), [4, 7, 0, 13, 8])


def test_construction_invalid_chars():
    with pytest.raises(Exception):
        DNAStringSet(["ACGT", "ACGT-X"])


def test_construction_empty():
    empty_set = DNAStringSet()
    assert len(empty_set) == 0
    assert empty_set.width().shape == (0,)

    empty_set_list = DNAStringSet([])
    assert len(empty_set_list) == 0
    assert empty_set_list.width().shape == (0,)


def test_names_setter(sample_set):
    new_names = ["a", "b", "c", "d", "e"]
    sample_set.names = new_names
    assert list(sample_set.names) == new_names

    with pytest.raises(Exception):
        sample_set.names = ["a", "b"]


def test_getitem_int(sample_set):
    seq = sample_set[1]
    assert isinstance(seq, DNAString)
    assert str(seq) == "GATTACA"

    empty_seq = sample_set[2]
    assert isinstance(empty_seq, DNAString)
    assert str(empty_seq) == ""

    last_seq = sample_set[-1]
    assert str(last_seq) == "ACGTACGT"


def test_getitem_slice_view(sample_set):
    subset = sample_set[1:4]
    assert isinstance(subset, DNAStringSet)
    assert len(subset) == 3
    assert list(subset.names) == ["seq2", "empty", "iupac"]
    assert np.array_equal(subset.width(), [7, 0, 13])

    # Test it's a view (shares the pool)
    assert id(subset._pool) == id(sample_set._pool)


def test_getitem_list_view(sample_set):
    indices = [0, 3, 4]
    subset = sample_set[indices]

    assert isinstance(subset, DNAStringSet)
    assert len(subset) == 3
    assert list(subset.names) == ["seq1", "iupac", "seq5"]
    assert np.array_equal(subset.width(), [4, 13, 8])

    # Test it's a view
    assert id(subset._pool) == id(sample_set._pool)


def test_to_list(sample_set, sample_seqs):
    assert sample_set.to_list() == [s.upper() for s in sample_seqs["seqs"]]


def test_unlist_contiguous(sample_set, sample_seqs):
    unlisted = sample_set.unlist()
    expected = "".join([s.upper() for s in sample_seqs["seqs"]])
    assert str(unlisted) == expected
    assert isinstance(unlisted, DNAString)


def test_unlist_subsetted_view(sample_set, sample_seqs):
    indices = [0, 2, 4]
    subset = sample_set[indices]

    unlisted = subset.unlist()
    expected = "ACGT" + "" + "ACGTACGT"
    assert str(unlisted) == expected
    assert isinstance(unlisted, DNAString)


def test_repr(sample_set):
    rep_str = repr(sample_set)
    assert "<DNAStringSet of length 5>" in rep_str
    assert "  [ 1]  4 ACGT                 seq1" in rep_str
    assert "  [ 3]  0                      empty" in rep_str
    assert "  [ 5]  8 ACGTACGT             seq5" in rep_str


def test_repr_long(sample_seqs):
    long_seqs = sample_seqs["seqs"] * 3
    long_names = sample_seqs["names"] * 3
    long_set = DNAStringSet(long_seqs, names=long_names)

    rep_str = repr(long_set)
    assert "<DNAStringSet of length 15>" in rep_str
    assert "..." in rep_str
    assert "  [ 1]  4 ACGT                 seq1" in rep_str
    assert "  [15]  8 ACGTACGT             seq5" in rep_str
    assert "  [ 5]  8 ACGTACGT             seq5" in rep_str
    assert "  [11]  4 ACGT                 seq1" in rep_str


def test_repr_empty():
    empty_set = DNAStringSet()
    assert repr(empty_set) == "<DNAStringSet of length 0>"
