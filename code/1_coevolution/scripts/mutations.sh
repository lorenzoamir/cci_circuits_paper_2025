#!/bin/bash
#PBS -N mutations
#PBS -l select=1:ncpus=16:ngpus=0:mem=16gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/1_coevolution
cd /home/lnemati/pathway_crosstalk/code/1_coevolution
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python get_mutations.py
exit 0
