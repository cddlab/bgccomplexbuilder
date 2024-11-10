import re
from dataclasses import dataclass
from pathlib import Path

from loguru import logger


@dataclass
class Predictedcomplex:
    name: str
    pLDDT: float
    pTM: float
    ipTM: float


def find_complexes_from_iptm(
    dirname: str | Path, threshold: float = 0.5
) -> list[Predictedcomplex]:
    """Find complexes from iptm score of ColabFold log.txt file.

    Args:
        - dirname (str | Path): Path to the directory containing ColabFold log.txt file.
        - threshold (float): Threshold for the confidence score. Default is 0.5.
    Returns:
        - complexes (list[Predictedcomplex]): List of predicted complexes.
    """
    dirname = Path(dirname)
    subdirs = [
        subdir.relative_to(dirname) for subdir in dirname.iterdir() if subdir.is_dir()
    ]
    # only rank_001
    # Example: BAV57458.1_BAV57459.1
    # 2024-11-09 11:15:18,557 rank_001_alphafold2_multimer_v3_model_4_seed_000
    # pLDDT=75.1 pTM=0.696 ipTM=0.215
    pattern = re.compile(
        r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3} "
        r"rank_001_alphafold2_multimer_v\d_model_\d+_seed_\d+ "
        r"pLDDT=(\d+\.\d) pTM=(\d+\.\d+) ipTM=(\d+\.\d+)"
    )

    results = []
    for subdir in subdirs:
        logfile = dirname / subdir / "log.txt"
        if not logfile.exists():
            logger.warning(f"{logfile} does not exist.")
            continue

        with open(logfile, "r") as f:
            for line in f:
                match = pattern.match(line)
                if match:
                    plddt, ptm, iptm = match.groups()
                    results.append(
                        Predictedcomplex(
                            name=str(subdir),
                            pLDDT=float(plddt),
                            pTM=float(ptm),
                            ipTM=float(iptm),
                        )
                    )
    threshold = 0.5
    complexes = [result for result in results if result.ipTM > threshold]

    return complexes
