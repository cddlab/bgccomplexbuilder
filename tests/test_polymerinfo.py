import pytest

from complexbuilder.polymerinfo.rcsbpdb_polymerinfo import (
    has_hetero_interfaces,
    has_homo_interfaces,
    make_assembly_interface_dict_from_pdbids,
)


@pytest.fixture(scope="module")
def assembly_rs():
    pdbids = ["3WWN", "7M7J", "13PK", "8QFU"]
    result_dict = make_assembly_interface_dict_from_pdbids(pdbids)
    return result_dict


def test_make_assembly_interface_dict_from_pdbids(assembly_rs):
    """Test make_assembly_interface_dict_from_pdbids with multiple cases."""
    assert (
        assembly_rs["entries"][0]["assemblies"][0]["interfaces"][0][
            "rcsb_interface_info"
        ]["interface_character"]
        == "homo"
    )


def test_has_homo_interfaces(assembly_rs):
    """Test has_homo_interfaces with multiple cases."""
    assert has_homo_interfaces(assembly_rs, "3WWN", "A") is True
    assert has_homo_interfaces(assembly_rs, "3WWN", "B") is False
    assert has_homo_interfaces(assembly_rs, "7M7J", "A") is True
    assert has_homo_interfaces(assembly_rs, "7M7J", "B") is True
    assert has_homo_interfaces(assembly_rs, "7M7J", "C") is False
    assert has_homo_interfaces(assembly_rs, "7M7J", "D") is False
    assert has_homo_interfaces(assembly_rs, "7M7J", "E") is False
    assert has_homo_interfaces(assembly_rs, "7M7J", "F") is False
    assert has_homo_interfaces(assembly_rs, "13PK", "A") is False
    assert has_homo_interfaces(assembly_rs, "13PK", "B") is False
    assert has_homo_interfaces(assembly_rs, "13PK", "C") is False
    assert has_homo_interfaces(assembly_rs, "13PK", "D") is False
    assert has_homo_interfaces(assembly_rs, "13PK", "E") is False  # Not exist
    assert has_homo_interfaces(assembly_rs, "8QFU", "A") is True
    assert has_homo_interfaces(assembly_rs, "8QFU", "B") is True
    assert has_homo_interfaces(assembly_rs, "8QFU", "C") is True
    assert has_homo_interfaces(assembly_rs, "8QFU", "D") is True
    assert has_homo_interfaces(assembly_rs, "8QFU", "E") is False  # Not exist


def test_has_hetero_interfaces(assembly_rs):
    """Test has_hetero_interfaces with multiple cases."""
    assert (
        has_hetero_interfaces(assembly_rs, "3WWN", "A", "A") is False
    )  # Homo dimer is not hetero
    assert has_hetero_interfaces(assembly_rs, "3WWN", "A", "B") is True
    assert (
        has_hetero_interfaces(assembly_rs, "3WWN", "A", "C") is False
    )  # chain C not exist
    assert has_hetero_interfaces(assembly_rs, "7M7J", "A", "B") is False  # Homo
    assert has_hetero_interfaces(assembly_rs, "7M7J", "B", "A") is False  # Homo
    assert has_hetero_interfaces(assembly_rs, "7M7J", "E", "F") is True
    assert has_hetero_interfaces(assembly_rs, "7M7J", "F", "E") is True
    assert has_hetero_interfaces(assembly_rs, "7M7J", "C", "D") is True
    assert has_hetero_interfaces(assembly_rs, "7M7J", "A", "E") is True
    assert has_hetero_interfaces(assembly_rs, "7M7J", "B", "C") is True
    assert has_hetero_interfaces(assembly_rs, "7M7J", "C", "D") is True
    assert has_hetero_interfaces(assembly_rs, "7M7J", "A", "D") is True
    assert has_hetero_interfaces(assembly_rs, "7M7J", "B", "F") is True
    assert has_hetero_interfaces(assembly_rs, "7M7J", "A", "F") is True
    assert has_hetero_interfaces(assembly_rs, "7M7J", "B", "D") is True
    assert has_hetero_interfaces(assembly_rs, "7M7J", "D", "E") is True
    assert has_hetero_interfaces(assembly_rs, "7M7J", "C", "F") is True
    assert has_hetero_interfaces(assembly_rs, "7M7J", "A", "H") is False  # Not exist
    assert (
        has_hetero_interfaces(assembly_rs, "7M7J", "A", "C") is True
    )  # But to be False
    assert has_hetero_interfaces(assembly_rs, "13PK", "A", "B") is False
    assert has_hetero_interfaces(assembly_rs, "8QFU", "A", "B") is False
    assert has_hetero_interfaces(assembly_rs, "8QFU", "C", "D") is False
