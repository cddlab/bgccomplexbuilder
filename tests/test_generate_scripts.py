import textwrap

import pytest

from complexbuilder.utils.generate_scripts import (
    generate_colabfold_search_runner,
)


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
    echo `hostname`

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
            "flow",
            "foo",
            "foo_dir",
            textwrap.dedent(
                """\
    #PJM -L rscunit=lm
    #PJM -L rscgrp=lm-middle
    #PJM -L socket=1
    #PJM -L elapse=2:00:00
    #PJM -j

    test $PJM_O_WORKDIR && cd $PJM_O_WORKDIR

    # run the environment module
    . /usr/share/Modules/init/sh
    module load gcc/11.3.0
    export PATH="/data/group1/z42347t/apps/localcolabfold/colabfold-conda/bin:$PATH"

    echo `hostname`

    colabfold_search \
      --use-env 1 \
      --use-templates 1 \
      --db-load-mode 2 \
      --db2 pdb100_230517 \
      --mmseqs /data/group1/z42347t/apps/mmseqs2/15-6f452/bin/mmseqs \
      --threads 4 \
      foo.fasta \
      /nvme1/z44243z/colabfold_env \
      foo_dir
    """
            ),
        ),
    ],
)
def test_generate_colabfold_search_runner(type, inputname, outputdir, expected_output):
    """test of generate_colabfold_search_runner"""
    actual_output = generate_colabfold_search_runner(type, inputname, outputdir)
    assert actual_output == expected_output
