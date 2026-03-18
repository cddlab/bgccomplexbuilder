# predicting protein complexes in biosynthetic gene clusters

This repository contains the code for a manuscript entitled "[Predicting Protein Complexes in Biosynthetic Gene Clusters](https://www.biorxiv.org/content/10.1101/2025.10.26.684697v1)".

The repository includes scripts for data preprocessing, evaluation, and visualization of protein complexes in biosynthetic gene clusters (BGCs). The resultant network maps and excel files generated from the analysis are available at Zenodo: [10.5281/zenodo.17451667](https://doi.org/10.5281/zenodo.17451667).

Please refer to the individual scripts for detailed instructions on how to use them.

## Installation

We employed [uv](https://docs.astral.sh/uv/) (developed by astral) and Python 3.12 for this project. First, install uv by following the instructions on the [official website](https://docs.astral.sh/uv/#installation). Then, clone this repository and set up a virtual environment:

```bash
git clone https://github.com/cddlab/complexbuilder.git
cd complexbuilder
uv sync
```

MiBIG GenBank files (GBK format) can be downloaded from the [MiBIG database](https://mibig.secondarymetabolites.org/download). We used MiBIG version 4.0 for this project. After downloading and extracting the GenBank files, you can use the `inputbuilder` script to create input FASTA files for `colabfold_search`.

## Usage

### Create input files for colabfold_search from GenBank files

```bash
uv run inputbuilder \
    -i ~/Downloads/mibig_gbk_4.0 \
    -o outputdir \
    --max_length 1950 \
    --maxbytes 6000000 \
    --start 1 \
    --end 2826
```

## Contact

For questions or further information, please contact the authors of the manuscript.
