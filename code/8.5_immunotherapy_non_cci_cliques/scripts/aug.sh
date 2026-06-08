#!/bin/bash
#PBS -N aug
#PBS -l select=1:ncpus=1:ngpus=0:mem=12gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/8_immunotherapy
cd /home/lnemati/pathway_crosstalk/code/8_immunotherapy
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate tabpfn
python augment.py
exit 0
