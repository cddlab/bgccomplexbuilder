# %%
import json
from pathlib import Path

import gemmi
from loguru import logger

from complexbuilder.common.log import log_setup

log_setup(level="INFO")


def calculate_rmsd_between_two_chains(
    ciffile: str | Path, chain1_id: str = "A", chain2_id: str = "B"
) -> float | None:
    """
    Compare two chains in a CIF file and return the RMSD value.

    Parameters:
        ciffile (str | Path): Path to the CIF file.
        chain1_id (str): ID of the first chain.
        chain2_id (str): ID of the second chain.

    Returns:
        float: RMSD value between the two chains.

    Note:
        If the chain length of either chain is less than 20, return "None"
        because the RMSD calculation is not reliable for short chains.
    """
    ciffile_str = str(ciffile) if isinstance(ciffile, Path) else ciffile
    if not Path(ciffile_str).exists():
        raise FileNotFoundError(f"CIF file not found: {ciffile_str}")
    model = gemmi.read_structure(ciffile_str)[0]
    chain1 = model[chain1_id].get_polymer()
    chain2 = model[chain2_id].get_polymer()

    # Ensure both chains are not empty
    assert chain1 is not None and chain2 is not None, "One of the chains is empty."
    # if the length of the chains is less than 20, return None
    if len(chain1) < 20 or len(chain2) < 20:
        return None

    ptype = chain1.check_polymer_type()
    sr = gemmi.calculate_superposition(
        chain1, chain2, ptype, gemmi.SupSelect.MainChain, trim_cycles=5
    )
    return round(sr.rmsd, 3)


def extract_hetero_complexes(json_data: dict) -> list[tuple[str, str]]:
    return [
        (bgc_id, complex_key)
        for bgc_id, complexes in json_data.items()
        for complex_key, complex_data in complexes.items()
        if complex_data.get("complex_homo_hetero") == "hetero"
    ]


# %%
bgcdirectory_root = "/data2/moriwaki/BGCcomplex/merged"
hitjsonfile = "/data2/moriwaki/complexbuilder/hitcomplex_iptm0.55_ipsae0.0_all.json"
rmsdoutput = "rmsd.json"

with open(hitjsonfile, "r") as f:
    hit_data = json.load(f)

count = 0
results = {}
with open(rmsdoutput, "w") as f:
    for bgc_id, complex_key in extract_hetero_complexes(hit_data):
        logger.info(f"Processing: {bgc_id}, {complex_key}")
        results[bgc_id] = results[bgc_id] if bgc_id in results else {}
        if Path(
            f"{bgcdirectory_root}/{bgc_id}/{complex_key}/{complex_key}_model.cif"
        ).exists():
            rmsd = calculate_rmsd_between_two_chains(
                f"{bgcdirectory_root}/{bgc_id}/{complex_key}/{complex_key}_model.cif"
            )
            if rmsd is not None:
                results[bgc_id][complex_key] = {"RMSD": rmsd}
        else:
            logger.warning(
                f"File not found: {bgcdirectory_root}/{bgc_id}/{complex_key}/{complex_key}_model.cif"
            )
            results[bgc_id][complex_key] = {"RMSD": None}

with open(rmsdoutput, "w") as f:
    json.dump(results, f, indent=2)
    logger.info(f"RMSD results saved to {rmsdoutput}")
# %%
