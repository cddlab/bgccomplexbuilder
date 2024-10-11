import pytest

from complexbuilder.utils.msa import make_multiple_msa


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
