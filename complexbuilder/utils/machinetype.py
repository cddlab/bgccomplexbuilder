import textwrap

yayoi = textwrap.dedent(
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
      {inputname}.fasta \
      /mnt/databases/colabfolddb \
      {outputdir}
    """
)

foodin = textwrap.dedent(
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
      {inputname}.fasta \
      /mnt/databases/colabfolddb \
      {outputdir}
      """
)

flow = textwrap.dedent(
    """\
    #!/bin/sh
    #PJM -L rscunit=lm
    #PJM -L rscgrp=lm-middle
    #PJM -L socket=1
    #PJM -L elapse=2:00:00
    #PJM -j

    test $PJM_O_WORKDIR && cd $PJM_O_WORKDIR

    # run the environment module
    . /usr/share/Modules/init/sh
    module load gcc/11.3.0
    export PATH="/data/group1/z44243z/apps/localcolabfold/colabfold-conda/bin:$PATH"

    echo `hostname`

    colabfold_search \
      --use-env 1 \
      --use-templates 1 \
      --db-load-mode 2 \
      --db2 pdb100_230517 \
      --mmseqs /data/group1/z44243z/mmseqs2/for_colabfold/bin/mmseqs \
      --threads 4 \
      {inputname}.fasta \
      /nvme1/z44243z/colabfold_env \
      {outputdir}
    """
)
