from complexbuilder.analysis.network import get_bgcgenes


def test_get_bgcgenes():
    dataInt = {
        "BGC0000001": {
            "gene1": {"name": "Gene One", "product": "Product One"},
            "gene2": {"name": "Gene Two", "product": "Product Two"},
        },
        "BGC0000002": {
            "gene3": {"name": "Gene Three", "product": "Product Three"},
            "gene4": {"name": "Gene Four", "product": "Product Four"},
        },
    }
    assert get_bgcgenes(dataInt, "BGC0000001") == ["gene1", "gene2"]
    assert get_bgcgenes(dataInt, "BGC0000002") == ["gene3", "gene4"]
