from complexbuilder.analysis.network import get_bgcgenes, make_node_desc_dict


def test_make_node_desc_dict():
    bgcgenes = ["aek75497.1", "aek75507.1", "abyU"]
    cds_info_list = [
        {
            "protein_id": "AEK75497.1",
            "gene": "gene1",
            "locus_tag": "locus1",
            "product": "product1",
        },
        {
            "protein_id": "AEK75507.1",
            "gene": None,
            "locus_tag": None,
            "product": "product2",
        },
        {
            "protein_id": None,
            "gene": "abyU",
            "locus_tag": None,
            "product": "product3",
        },
    ]
    expected_result = {
        "aek75497.1": "AEK75497.1\ngene1\nlocus1\nproduct1",
        "aek75507.1": "AEK75507.1\nproduct2",
        "abyU": "abyU\nproduct3",
    }
    assert make_node_desc_dict(bgcgenes, cds_info_list) == expected_result
