import textwrap

import pytest

from complexbuilder.utils.generate_scripts import generate_scripts


@pytest.mark.parametrize(
    "type, inputname, outputdir, expected_output",
    [
        pytest.param(
            "yayoi",
            "foo",
            "foo_dir",
            textwrap.dedent(
                """\
    #!/bin/sh
    #PBS -l nodes=1:ppn=4:yayoims
    #PBS -l walltime=72:00:00
    #PBS -q default

    test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR

    # run the environment module
    . /home/apps/Modules/init/profile.sh
    module load localcolabfold
    colabfold_search \
      --use-env 1 \
      --use-templates 1 \
      --db-load-mode 2 \
      --db2 pdb100_230517 \
      --mmseqs /home/apps/mmseqs2/15-6f452/bin/mmseqs \
      --threads 4 \
      foo.fasta \
      /mnt/databases/colabfolddb \
      foo_dir
    """
            ),
        ),
        pytest.param(
            "foodin",
            "bar",
            "bar_dir",
            textwrap.dedent(
                """\
    #!/bin/sh
    #SBATCH -p all_q
    #SBATCH -n 32
    #SBATCH --gpus 1
    #SBATCH --partition=hf
    #SBATCH --time=04:00:00
    #SBATCH -o %x.%j.out
    #SBATCH -e %x.%j.err

    . /home/apps/Modules/init/profile.sh
    module load cuda/11.8 gcc/13.3.0

    colabfold_search \
      --use-env 1 \
      --use-templates 1 \
      --db-load-mode 2 \
      --db2 pdb100_230517 \
      --mmseqs /home/apps/mmseqs2/15-6f452/bin/mmseqs \
      --threads 4 \
      bar.fasta \
      /mnt/databases/colabfolddb \
      bar_dir
      """
            ),
        ),
    ],
)
def test_generate_scripts(type, inputname, outputdir, expected_output):
    """test of generate_scripts"""
    assert generate_scripts(type, inputname, outputdir) == expected_output
