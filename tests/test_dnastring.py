import pytest

from biostrings import DNAString


def test_construction_str():
    seq_str = "ACGT-N"
    dna_str = DNAString(seq_str)
    assert len(dna_str) == 6
    assert str(dna_str) == "ACGT-N"
    assert dna_str.to_bytes() == b"ACGT-N"


def test_construction_bytes():
    seq_bytes = b"ACGT-N"
    dna_str = DNAString(seq_bytes)
    assert len(dna_str) == 6
    assert str(dna_str) == "ACGT-N"
    assert dna_str.to_bytes() == b"ACGT-N"


def test_construction_lowercase_and_validation():
    dna_str = DNAString("acgt-n")
    assert str(dna_str) == "ACGT-N"


def test_construction_invalid_chars():
    with pytest.raises(ValueError, match="non-DNA characters"):
        DNAString("ACGT-X")


def test_construction_invalid_type():
    with pytest.raises(TypeError):
        DNAString(123)


def test_len():
    assert len(DNAString("ACGT")) == 4
    assert len(DNAString("")) == 0


def test_repr():
    dna = DNAString("ACGT")
    assert repr(dna) == "DNAString(length=4, sequence='ACGT')"
    long_dna = DNAString("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    assert "..." in repr(long_dna)
    assert repr(long_dna) == "DNAString(length=29, sequence='AAAAAAAAAA...AAAAAAAAAA')"


def test_eq():
    dna1 = DNAString("ACGT")
    dna2 = DNAString("ACGT")
    dna3 = DNAString("TCGA")
    assert dna1 == dna2
    assert dna1 != dna3
    assert dna1 == "ACGT"
    assert dna1 == "acgt"
    assert dna1 != "TCGA"
    assert dna1 != 123


def test_getitem_int():
    dna = DNAString("ACGT")
    assert dna[0] == "A"
    assert dna[1] == "C"
    assert dna[-1] == "T"
    assert isinstance(dna[0], DNAString)


def test_getitem_slice():
    dna = DNAString("AAACCCGGG")
    sub = dna[3:6]
    assert isinstance(sub, DNAString)
    assert str(sub) == "CCC"

    sub_all = dna[:]
    assert str(sub_all) == "AAACCCGGG"


def test_reverse_complement_simple():
    dna = DNAString("AAACCCGGGTTT")
    rc = dna.reverse_complement()
    assert str(rc) == "AAACCCGGGTTT"


def test_reverse_complement_iupac():
    #
    dna = DNAString("ACGTRYSWKMBDHVN-")
    rc = dna.reverse_complement()
    # From our complement table
    assert str(rc) == "-NBDHVKMWSRYACGT"
