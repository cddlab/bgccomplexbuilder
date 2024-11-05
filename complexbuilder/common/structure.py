from pathlib import Path


def truncate_low_plddt_terminal_residues(pdbfile: Path, threshold: float = 50.0):
    """Truncate terminal residues with low PLDDT scores.
    inputs:
    - pdbfile: Path to the PDB file.
    - threshold: Threshold for truncating terminal residues.
    returns:
    - structure: Bio.PDB.Structure, a structure object
    """
    # assert pdbfile.exists()
    if not pdbfile.exists():
        raise FileNotFoundError(f"{pdbfile} does not exist.")
