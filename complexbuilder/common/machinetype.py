import dataclasses
import textwrap


@dataclasses.dataclass
class yayoi:
    dbpath: str = "/mnt/databases/colabfolddb"
    localpdbpath: str = "/home/database/pdb_mmcif/mmcif_files"

    @classmethod
    def colabfold_search(cls, inputname: str, outputdir: str) -> str:
        return textwrap.dedent(
            f"""\
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
      {inputname}.fasta \
      {cls.dbpath} \
      {outputdir}
    """
        )

    @classmethod
    def colabfold_batch(cls, inputname: str, outputdir: str) -> str:
        return textwrap.dedent(
            f"""\
    #!/bin/sh
    #PBS -q default
    #PBS -l nodes=1:ppn=16:gpus=1
    #PBS -l walltime=24:00:00

    test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR

    # run the environment module
    . /home/apps/Modules/init/profile.sh
    module load localcolabfold
    echo `hostname`

    RANDOMSEED=0
    OUTPUTDIR={outputdir}
    INPUTFILE="{inputname}.a3m"
    PDBHITFILE="{inputname}_pdb100_230517.m8"
    LOCALPDBPATH="{cls.localpdbpath}"

    # set recycle to 1 to perform only one recycle
    colabfold_batch \
      --num-recycle 1 \
      --amber \
      --templates \
      --use-gpu-relax \
      --num-models 2 \
      --model-order 1,2 \
      --pdb-hit-file ${{PDBHITFILE}} \
      --local-pdb-path ${{LOCALPDBPATH}} \
      --random-seed ${{RANDOMSEED}} \
      ${{INPUTFILE}} \
      ${{OUTPUTDIR}}
      """
        )


class flow:
    dbpath: str = "/nvme1/z44243z/colabfold_env"
    localpdbpath: str = "/beegfs/share/alphafold/db-v2.3/pdb_mmcif/mmcif_files"

    @classmethod
    def colabfold_search(cls, inputname: str, outputdir: str) -> str:
        return textwrap.dedent(
            f"""\
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
      {inputname}.fasta \
      {cls.dbpath} \
      {outputdir}
    """
        )

    @classmethod
    def colabfold_batch(cls, inputname: str, outputdir: str) -> str:
        return textwrap.dedent(
            f"""\
    #!/bin/sh
    #PJM -L rscunit=cx
    #PJM -L rscgrp=cxgfs-single
    #PJM -L node=1
    #PJM -L elapse=01:00:00
    #PJM -j
    #PJM -S

    test $PJM_O_WORKDIR && cd $PJM_O_WORKDIR

    # run the environment module
    . /usr/share/Modules/init/sh
    module use --append /data/group1/z42347t/modulefiles
    module load localcolabfold

    echo `hostname`

    RANDOMSEED=0
    OUTPUTDIR={outputdir}
    INPUTFILE="{inputname}.a3m"
    PDBHITFILE="{inputname}_pdb100_230517.m8"
    LOCALPDBPATH="{cls.localpdbpath}"

    # set recycle to 1 to perform only one recycle
    colabfold_batch \
      --num-recycle 1 \
      --amber \
      --templates \
      --use-gpu-relax \
      --num-models 2 \
      --model-order 1,2 \
      --pdb-hit-file ${{PDBHITFILE}} \
      --local-pdb-path ${{LOCALPDBPATH}} \
      --random-seed ${{RANDOMSEED}} \
      ${{INPUTFILE}} \
      ${{OUTPUTDIR}}
      """
        )
