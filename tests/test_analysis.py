import json
from pathlib import Path

# from complexbuilder.analysis.conjoinedtwins import (
#     calculate_rmsd_between_two_chains,
#     extract_hetero_complexes,
# )
from complexbuilder.analysis.extract_hitcomplexes import (
    _check_homo_hetero,
    make_hitcomplexlist,
    split_proteinids,
)

# def test_calculate_rmsd_between_two_chains():
#     ciffile = Path("tests/testfiles/MonBI_MonBII_model.cif")
#     rmsd = calculate_rmsd_between_two_chains(ciffile, "A", "B")
#     assert rmsd is not None and rmsd < 0.7


# def test_extract_hetero_complexes():
#     test_data = {
#         "BGC0000001": {
#             "aek75492.1_aek75495.1": {
#                 "complex_homo_hetero": "hetero",
#             }
#         },
#         "BGC0000007": {
#             "aas90001.1_aas90002.1": {
#                 "complex_homo_hetero": "hetero",
#             },
#             "aas90020.1_aas90020.1": {
#                 "complex_homo_hetero": "homo",
#             },
#             "aas89996.1_aas89996.1": {
#                 "complex_homo_hetero": "homo",
#             },
#             "aas90003.1_aas90004.1": {
#                 "complex_homo_hetero": "hetero",
#             },
#         },
#     }
#     assert extract_hetero_complexes(test_data) == [
#         ("BGC0000001", "aek75492.1_aek75495.1"),
#         ("BGC0000007", "aas90001.1_aas90002.1"),
#         ("BGC0000007", "aas90003.1_aas90004.1"),
#     ]


def test_check_homo_hetero():
    assert _check_homo_hetero("adakdk.aa_adakdk.aa") == "homo"
    assert _check_homo_hetero("alsk.01_alsk.01") == "homo"
    assert _check_homo_hetero("akskk_002.01_akskk_002.01") == "homo"
    assert _check_homo_hetero("allsk_01.002_allsk_01.003") == "hetero"
    assert _check_homo_hetero("trx17522.1_trx20192.1") == "hetero"


def test_split_proteinids():
    assert split_proteinids("trx17522.1_trx20192.1") == ("trx17522.1", "trx20192.1")
    # include "xx_04290"
    assert split_proteinids("fnf07_04290_trx17522.1") == (
        "fnf07_04290",
        "trx17522.1",
    )
    assert split_proteinids("trx17522.1_fnf07_04290") == (
        "trx17522.1",
        "fnf07_04290",
    )
    # case of "_orf"
    assert split_proteinids("ctg1_orf2_ctg1_orf2") == (
        "ctg1_orf2",
        "ctg1_orf2",
    )
    assert split_proteinids("wp_032798144.1_ssgg_rs34700") == (
        "wp_032798144.1",
        "ssgg_rs34700",
    )
    assert split_proteinids("rso11565.1_rso11564.1") == (
        "rso11565.1",
        "rso11564.1",
    )
    assert split_proteinids("ralta_b1233_caq71835.1") == (
        "ralta_b1233",
        "caq71835.1",
    )
    assert split_proteinids("ralta_b1233_ralta_b1233") == (
        "ralta_b1233",
        "ralta_b1233",
    )
    assert split_proteinids("aerm__aerm_") == ("aerm_", "aerm_")
    assert split_proteinids("acm68692.1_aerm_") == ("acm68692.1", "aerm_")
    assert split_proteinids("aerm__acm68692.1") == ("aerm_", "acm68692.1")


def test_make_hitcomplexlist(tmp_path):
    # Setup directory for BGC0000001
    dir1 = tmp_path / "BGC0000001"
    dir1.mkdir()
    metrics1 = {
        "ipSAE": [
            {"aek75497.1_aek75497.1": 0.918005},
            {"aek75507.1_aek75507.1": 0.87621},
            {"aek75514.1_aek75514.1": 0.864171},
            {"abyu_abyu": 0.84221},
            {"aek75504.1_aek75504.1": 0.746544},
            {"aek75492.1_aek75495.1": 0.742241},
        ],
        "ipTM": [
            {"aek75497.1_aek75497.1": 0.95},
            {"aek75507.1_aek75507.1": 0.92},
            {"aek75514.1_aek75514.1": 0.89},
            {"abyu_abyu": 0.89},
            {"aek75504.1_aek75504.1": 0.8},
            {"aek75492.1_aek75495.1": 0.85},
        ],
    }
    with open(dir1 / "complexmetrics.json", "w") as f:
        json.dump(metrics1, f)

    # Setup directory for BGC0000002
    dir2 = tmp_path / "BGC0000002"
    dir2.mkdir()
    metrics2 = {
        "ipSAE": [
            {"ahh99946.1_ahh99946.1": 0.892961},
            {"ahh99942.1_ahh99942.1": 0.879215},
            {"ahh99937.1_ahh99937.1": 0.853516},
            {"ahh99935.1_ahh99935.1": 0.789216},
            {"ahh99932.1_ahh99932.1": 0.780463},
            {"ahh99944.1_ahh99944.1": 0.774838},
            {"ahh99940.1_ahh99940.1": 0.76713},
            {"ahh99941.1_ahh99941.1": 0.758867},
            {"ahh99939.1_ahh99939.1": 0.715651},
            {"ahh99930.1_ahh99930.1": 0.662392},
            {"ahh99918.1_ahh99936.1": 0.653957},
        ],
        "ipTM": [
            {"ahh99946.1_ahh99946.1": 0.93},
            {"ahh99942.1_ahh99942.1": 0.93},
            {"ahh99937.1_ahh99937.1": 0.89},
            {"ahh99935.1_ahh99935.1": 0.86},
            {"ahh99932.1_ahh99932.1": 0.83},
            {"ahh99944.1_ahh99944.1": 0.87},
            {"ahh99940.1_ahh99940.1": 0.88},
            {"ahh99941.1_ahh99941.1": 0.85},
            {"ahh99939.1_ahh99939.1": 0.85},
            {"ahh99930.1_ahh99930.1": 0.84},
            {"ahh99918.1_ahh99936.1": 0.81},
        ],
    }
    with open(dir2 / "complexmetrics.json", "w") as f:
        json.dump(metrics2, f)

    # Create extra directory without complexmetrics.json to test skipping.
    extra_dir = tmp_path / "BGC0000003"
    extra_dir.mkdir()

    # Run the function
    result = json.dumps(make_hitcomplexlist(tmp_path, 0.6, 0.8), indent=4)
    result_dict = json.loads(result)
    expected = {
        "BGC0000001": {
            "aek75497.1_aek75497.1": {
                "ipSAE": 0.918005,
                "ipTM": 0.95,
                "complex_homo_hetero": "homo",
            },
            "aek75507.1_aek75507.1": {
                "ipSAE": 0.87621,
                "ipTM": 0.92,
                "complex_homo_hetero": "homo",
            },
            "aek75514.1_aek75514.1": {
                "ipSAE": 0.864171,
                "ipTM": 0.89,
                "complex_homo_hetero": "homo",
            },
            "abyu_abyu": {
                "ipSAE": 0.84221,
                "ipTM": 0.89,
                "complex_homo_hetero": "homo",
            },
            "aek75504.1_aek75504.1": {
                "ipSAE": 0.746544,
                "ipTM": 0.8,
                "complex_homo_hetero": "homo",
            },
            "aek75492.1_aek75495.1": {
                "ipSAE": 0.742241,
                "ipTM": 0.85,
                "complex_homo_hetero": "hetero",
            },
        },
        "BGC0000002": {
            "ahh99946.1_ahh99946.1": {
                "ipSAE": 0.892961,
                "ipTM": 0.93,
                "complex_homo_hetero": "homo",
            },
            "ahh99942.1_ahh99942.1": {
                "ipSAE": 0.879215,
                "ipTM": 0.93,
                "complex_homo_hetero": "homo",
            },
            "ahh99937.1_ahh99937.1": {
                "ipSAE": 0.853516,
                "ipTM": 0.89,
                "complex_homo_hetero": "homo",
            },
            "ahh99935.1_ahh99935.1": {
                "ipSAE": 0.789216,
                "ipTM": 0.86,
                "complex_homo_hetero": "homo",
            },
            "ahh99932.1_ahh99932.1": {
                "ipSAE": 0.780463,
                "ipTM": 0.83,
                "complex_homo_hetero": "homo",
            },
            "ahh99944.1_ahh99944.1": {
                "ipSAE": 0.774838,
                "ipTM": 0.87,
                "complex_homo_hetero": "homo",
            },
            "ahh99940.1_ahh99940.1": {
                "ipSAE": 0.76713,
                "ipTM": 0.88,
                "complex_homo_hetero": "homo",
            },
            "ahh99941.1_ahh99941.1": {
                "ipSAE": 0.758867,
                "ipTM": 0.85,
                "complex_homo_hetero": "homo",
            },
            "ahh99939.1_ahh99939.1": {
                "ipSAE": 0.715651,
                "ipTM": 0.85,
                "complex_homo_hetero": "homo",
            },
            "ahh99930.1_ahh99930.1": {
                "ipSAE": 0.662392,
                "ipTM": 0.84,
                "complex_homo_hetero": "homo",
            },
            "ahh99918.1_ahh99936.1": {
                "ipSAE": 0.653957,
                "ipTM": 0.81,
                "complex_homo_hetero": "hetero",
            },
        },
    }
    assert result_dict == expected
