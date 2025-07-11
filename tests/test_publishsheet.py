import pandas as pd

from complexbuilder.pdbinbgcclusters import homomer_add_hit_to_df


def test_homomer_add_hit_to_df():
    """Test the homomer_add_hit_to_df function with sample data."""
    # Sample data
    sample_data = {
        "BGC0000017": {"acr33074.1_4mer": {"ipSAE": 0.48, "ipTM": 0.40}},
        "BGC0000033": {
            "aam94792.1_4mer": {"ipSAE": 0.57, "ipTM": 0.33},
            "aam94793.1_4mer": {"ipSAE": 0.83, "ipTM": 0.19},
            "aam70337.1_10mer": {"ipSAE": 0.116, "ipTM": 0.893},
        },
        "BGC0000081": {},
        "BGC0000082": {"wp_003722052.1_3mer": {"ipSAE": 0.011, "ipTM": 0.89}},
    }

    df = homomer_add_hit_to_df(sample_data)
    assert isinstance(df, pd.DataFrame), "Output should be a DataFrame"
    assert df["mibig_accession"][0] == "BGC0000017"
    assert df["proteins"][0] == "acr33074.1_acr33074.1"
    assert df["ipSAE"][0] == 0.48
    assert df["ipTM"][0] == 0.40
    assert df["mibig_accession"][1] == "BGC0000033"
    assert df["proteins"][1] == "aam94792.1_aam94792.1"
    assert df["ipSAE"][1] == 0.57
    assert df["ipTM"][1] == 0.33
    assert df["mibig_accession"][2] == "BGC0000033"
    assert df["proteins"][2] == "aam94793.1_aam94793.1"
    assert df["ipSAE"][2] == 0.83
    assert df["ipTM"][2] == 0.19
    assert df["mibig_accession"][3] == "BGC0000033"
    assert df["proteins"][3] == "aam70337.1_aam70337.1"
    assert df["ipSAE"][3] == 0.116
    assert df["ipTM"][3] == 0.893
    assert df["mibig_accession"][4] == "BGC0000082"
    assert df["proteins"][4] == "wp_003722052.1_wp_003722052.1"
    assert df["ipSAE"][4] == 0.011
    assert df["ipTM"][4] == 0.89
    assert len(df) == 5, "DataFrame should have 5 rows"
