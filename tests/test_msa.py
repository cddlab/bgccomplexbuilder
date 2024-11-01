import pytest

from complexbuilder.utils.msa import (
    get_protein_sequence_from_uniprot,
    make_multiple_msa,
)


@pytest.mark.parametrize(
    "seqs, copies, expected_output",
    [
        pytest.param(
            {"seq1": "MTEITAAMVKELREST", "seq2": "AKAIKES"},
            (1, 2),
            "MTEITAAMVKELREST:AKAIKES:AKAIKES",
            id="positive case",
        ),
        pytest.param(
            {"seq3": "MTEITAAMVKELRESTAA"},
            (2,),
            "MTEITAAMVKELRESTAA:MTEITAAMVKELRESTAA",
            id="test2",
        ),
        pytest.param(
            {"seq3": "MTEITAAMVKELRESTAA"},
            2,
            "MTEITAAMVKELRESTAA:MTEITAAMVKELRESTAA",
            id="test3",
        ),
        pytest.param(
            {"seq1": "MTEITAAMVKELREST", "seq2": "AKAIKES"},
            (2, 1),
            "MTEITAAMVKELREST:MTEITAAMVKELREST:AKAIKES",
            id="test4",
        ),
    ],
)
def test_make_multiple_msa(seqs, copies, expected_output):
    """test of make_multiple_msa"""
    assert make_multiple_msa(seqs, copies).seq == expected_output
    assert make_multiple_msa(seqs, copies).id == "_".join(seqs.keys())


@pytest.mark.parametrize(
    "uniprot_id, expected_output, raises_exception",
    [
        (
            "P69905",
            (
                "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPN"
                "ALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
            ),
            False,  # 正常系: 例外が発生しない
        ),
        (
            "P9C9999",
            None,  # No expected output
            True,  # Invalid uniprot_id
        ),
    ],
)
def test_get_protein_sequence_from_uniprot(
    uniprot_id, expected_output, raises_exception
):
    """Test get_protein_sequence_from_uniprot with multiple cases."""
    if raises_exception:
        with pytest.raises(ValueError, match="Failed to retrieve data"):
            get_protein_sequence_from_uniprot(uniprot_id)
    else:
        sequence = get_protein_sequence_from_uniprot(uniprot_id)
        assert (
            sequence == expected_output
        ), "The sequence does not match the expected result."
