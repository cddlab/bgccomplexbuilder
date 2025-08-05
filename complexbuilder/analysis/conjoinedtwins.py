# %%
from pathlib import Path

import gemmi


def compare_two_chains(
    ciffile: str | Path, chain1_id: str = "A", chain2_id: str = "B"
) -> float:
    """
    Compare two chains in a CIF file and return the RMSD value.

    Parameters:
        ciffile (str | Path): Path to the CIF file.
        chain1_id (str): ID of the first chain.
        chain2_id (str): ID of the second chain.

    Returns:
        float: RMSD value between the two chains.
    """
    ciffile_str = str(ciffile) if isinstance(ciffile, Path) else ciffile
    if not Path(ciffile_str).exists():
        raise FileNotFoundError(f"CIF file not found: {ciffile_str}")
    model = gemmi.read_structure(ciffile_str)[0]
    chain1 = model[chain1_id].get_polymer()
    chain2 = model[chain2_id].get_polymer()

    # Ensure both chains are not empty
    assert chain1 is not None and chain2 is not None, "One of the chains is empty."

    ptype = chain1.check_polymer_type()
    sr = gemmi.calculate_superposition(
        chain1, chain2, ptype, gemmi.SupSelect.MainChain, trim_cycles=5
    )
    return sr.rmsd


# %%
