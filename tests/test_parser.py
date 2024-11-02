import pytest

from complexbuilder.utils.parser import classify_proteins, parse_mibig_json


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


def test_classify_proteins():
    """Test classify_proteins with multiple cases."""
    genbank_file = "/Users/YoshitakaM/Downloads/mibig_gbk_3.1/BGC0000028.gbk"
    nonnrpspksproteins, nrpspksproteins = classify_proteins(genbank_file)
    assert nonnrpspksproteins[0].id == "ADC79613.1"
    expected_seq = "MTLSVASILSESALRRPEHPAVVSGTRKTTYRELWDEARRYAAALRARGIGPGDKVALLLPSTPHFPSAYFGVLALGAIAVPVHALLRADEIAYILKDSGAAALICAAPLLAEGGRAAETTGTPVFTVMAERDEAARASAPRLDALAARSTPIDRQVPRAPEDIAVILYTSGTTGRPKGALLTHLNVVMNVDTTMLSPFDFTADDVLLGCLPLFHTFGQICGMNTCFRAGATLVLMPRFDGPDALDLLVREGCTVFMGVPTMYTALLEAARADPRRPALDRAFSGGAALPVAVLDAFRETFGCPVLEGYGLTETSPVVAYNQRAWPLRPGTVGRPIWGVEVEIARAEVEDRIELLPVGETGEIVIRGHNVMAGYLNRPEATAEAIVDGWFRSGDLGVKDDEGYLSVVDRKKDVVLRGGYNVYPREVEDVLAHHPAIAQAAVVGLPHPVHGEEVCAVVRPHPGTAPDPALGAEIVAWSKERMAPYKYPRRVEFVDAFPLGPSGKVLKRELVARLTAGARQVRTEAQETA"  # noqa
    assert nonnrpspksproteins[0].seq == expected_seq
    assert nonnrpspksproteins[0].description == "BafX"
    assert len(nonnrpspksproteins) == 13
    assert nrpspksproteins[0].id == "ADC79616.1"
    assert nrpspksproteins[0].description == "ADC79616.1_PKS_AT.1"
    assert len(nrpspksproteins[0].seq) == 299
    expected_seq = "VFPGQGAQWPRMAVDLLDTSTVFRDRMDACAQALEPFVDWSPLDVLRAPGAPGAPAGDRADVVQPLLFAVTVSLAALWRSHGVEPAAVLGHSVGEVTAAAVSGALSLDDSARVVALWSQAQATLAGQGDMVSVMAPAAEVEPRLHRWEGRLVVAAHNGPRSVIVSGDRDAAAELLDGLAADAVHARRIAVGLAAHSPHIDAIIPRMRADLAPIHPRTPHLPYYSGLTGGRLDAPALDADYWCRNLRNTVRFHQAARALLRDGHGVLLEVSPHTVLTSALTDCVEEHGVQAAVLGTLRRD"  # noqa
    assert nrpspksproteins[0].seq == expected_seq
    assert len(nrpspksproteins) == 59
