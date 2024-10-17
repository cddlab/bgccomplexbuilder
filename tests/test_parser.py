from complexbuilder.utils.parser import parse_mibig_json


def test_parse_mibig_json():
    """Test get_protein_sequence with multiple cases."""
    mibig_json = "/Users/YoshitakaM/Downloads/mibig_json_3.1/BGC0000028.json"
    mibig_data = parse_mibig_json(mibig_json)
    assert mibig_data["cluster"]["mibig_accession"] == "BGC0000028"
    assert mibig_data["cluster"]["biosyn_class"] == ["Polyketide"]
    assert mibig_data["cluster"]["compounds"][0]["compound"] == "bafilomycin B1"
    assert mibig_data["cluster"]["compounds"][0]["chem_struct"].startswith(
        "COC1\\C=C\\C=C(C)"
    )
