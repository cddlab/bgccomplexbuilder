import Bio.PDB.alphafold_db
from pathlib import Path


def truncate_low_plddt_terminal_residues(pdbfile: Path, threshold=0.5):
    """Truncate terminal residues with low PLDDT scores."""
    # assert pdbfile.exists()
    if not pdbfile.exists():
        raise FileNotFoundError(f"{pdbfile} does not exist.")

    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdbfile)
