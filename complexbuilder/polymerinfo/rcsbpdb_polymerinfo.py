# %%
from loguru import logger
from rcsbapi.data import DataQuery as Query

from complexbuilder.common.log import log_setup

log_setup(level="INFO")


def make_assembly_interface_dict_from_pdbids(pdbids: list[str]) -> dict:
    """Make a dict of interfaces from a list of PDB IDs.
    This function queries the RCSB PDB API to retrieve assembly interface information
    for a list of PDB IDs. It returns a dictionary where each key is a PDB ID and the value
    is a list of interfaces associated with that PDB ID.

    Args:
        pdbids (list): List of PDB IDs. e.g. ["3WWN", "7M7J", "13PK", "8QFU"]
    Returns:
        dict: A dictionary with PDB IDs as keys and a list of interfaces as values.
    """
    query = Query(
        input_type="entries",
        input_ids=pdbids,
        return_data_list=[
            "rcsb_id",
            "rcsb_entry_container_identifiers.polymer_entity_ids",
            "polymer_entities.entity_poly.pdbx_strand_id",
            "assemblies.interfaces.rcsb_interface_partner.interface_partner_identifier",
            "assemblies.interfaces.rcsb_interface_info.interface_character",
        ],
    )
    return query.exec()["data"]


# %%
def has_homo_interfaces(interface_dict: dict, pdbid: str, chainid: str) -> bool:
    """
    Check if a given PDB ID and chain ID has homo interfaces in the provided interface dictionary.

    This function searches for the specified PDB ID within the interface dictionary, verifies whether
    the given chain ID is present in the corresponding polymer entities (by comparing against the
    pdbx_strand_id list), and checks if any interface associated with the determined entity number is
    characterized as "homo".

    Args:
        interface_dict (dict): A dictionary containing data from a query to the RCSB PDB API. Expected
                               to have an "entries" key with a list of entry dictionaries.
        pdbid (str): The PDB ID to search for in the interface dictionary.
        chainid (str): The chain identifier to check in the polymer entities for the specified PDB ID.

    Returns:
        bool: True if a homo interface is found for the specified PDB ID and chain ID; False otherwise.
    """
    # Check if 'pdbid' exists in interface_dict["entries"][*]["rcsb_id"]
    entry_no = None
    for i, entry in enumerate(interface_dict["entries"]):
        if entry["rcsb_id"] == pdbid:
            entry_no = i
            break

    if entry_no is None:
        logger.debug(f"PDB ID: {pdbid} not found in interface dictionary.")
        return False
    else:
        logger.debug(f"Found PDB ID: {pdbid} in interface dictionary.")
        # If "pdbx_strand_id" is "B,A" for example, split it by "," and check if it contains chainid (i.e. ["B", "A"])
        entity_no = None
        for j, polymer_entity in enumerate(
            interface_dict["entries"][entry_no]["polymer_entities"]
        ):
            if chainid in polymer_entity["entity_poly"]["pdbx_strand_id"].split(","):
                entity_no = j + 1
                break
        if entity_no is None:
            logger.debug(
                f"Chain ID: {chainid} not found in PDB entry: {pdbid} in interface dictionary."
            )
            return False
        else:
            logger.debug(
                f"Found Chain ID: {chainid} in PDB ID: {pdbid} at entity ID: {entity_no}."
            )
            # Check if entity_no exists in the interfaces of the assemblies
            # interface_dict["entries"][entry_no]["assemblies"][*]["interfaces"][*]["rcsb_interface_partner"][*]["interface_partner_identifier"]["entity_id"]
            for assembly in interface_dict["entries"][entry_no]["assemblies"]:
                if assembly["interfaces"] is None:
                    logger.debug(
                        f"No interfaces found in assembly for PDB ID: {pdbid}."
                    )
                    return False
                logger.debug(f"Checking interfaces in assembly for PDB ID: {pdbid}.")
                for interface in assembly["interfaces"]:
                    for partner in interface["rcsb_interface_partner"]:
                        if (
                            partner["interface_partner_identifier"]["entity_id"]
                            == str(entity_no)
                            and interface["rcsb_interface_info"]["interface_character"]
                            == "homo"
                        ):
                            return True
    return False


# %%
def has_hetero_interfaces(
    interface_dict: dict, pdbid: str, chainid_1: str, chainid_2: str
) -> bool:
    """
    Check if a given PDB ID and chain ID has hetero interfaces in the provided interface dictionary.
    """
    entry_no = None
    for i, entry in enumerate(interface_dict["entries"]):
        if entry["rcsb_id"] == pdbid:
            entry_no = i
            break

    if entry_no is None:
        logger.debug(f"PDB ID: {pdbid} not found in interface dictionary.")
    else:
        chainid_1_entity_no = None
        chainid_2_entity_no = None
        for j, polymer_entity in enumerate(
            interface_dict["entries"][entry_no]["polymer_entities"]
        ):
            if chainid_1 in polymer_entity["entity_poly"]["pdbx_strand_id"].split(","):
                chainid_1_entity_no = j + 1
            if chainid_2 in polymer_entity["entity_poly"]["pdbx_strand_id"].split(","):
                chainid_2_entity_no = j + 1
        if chainid_1_entity_no is None or chainid_2_entity_no is None:
            logger.debug(
                f"Chain ID: {chainid_1} or {chainid_2} not found in PDB entry: {pdbid} in interface dictionary."
            )
            return False
        else:
            logger.debug(
                f"Found Chain ID: {chainid_1} in PDB ID: {pdbid} at entity ID: {chainid_1_entity_no}."
            )
            # Check if entity_no exists in the interfaces of the assemblies
            # interface_dict["entries"][entry_no]["assemblies"][*]["interfaces"][*]["rcsb_interface_partner"][{0,1}]["interface_partner_identifier"]["entity_id"]
            for assembly in interface_dict["entries"][entry_no]["assemblies"]:
                if assembly["interfaces"] is None:
                    logger.debug(
                        f"No interfaces found in assembly for PDB ID: {pdbid}."
                    )
                    return False
                logger.debug(f"Checking interfaces in assembly for PDB ID: {pdbid}.")
                for interface in assembly["interfaces"]:
                    logger.debug(f"Checking interface: {interface}.")
                    if (
                        interface["rcsb_interface_partner"][0][
                            "interface_partner_identifier"
                        ]["entity_id"]
                        == str(chainid_1_entity_no)
                        and interface["rcsb_interface_partner"][1][
                            "interface_partner_identifier"
                        ]["entity_id"]
                        == str(chainid_2_entity_no)
                    ) or (
                        interface["rcsb_interface_partner"][0][
                            "interface_partner_identifier"
                        ]["entity_id"]
                        == str(chainid_2_entity_no)
                        and interface["rcsb_interface_partner"][1][
                            "interface_partner_identifier"
                        ]["entity_id"]
                        == str(chainid_1_entity_no)
                    ):
                        if (
                            interface["rcsb_interface_info"]["interface_character"]
                            == "hetero"
                        ):
                            logger.debug(
                                f"Hetero interface found for PDB ID: {pdbid}, Chain IDs: {chainid_1}, {chainid_2}."
                            )
                            return True
    return False
