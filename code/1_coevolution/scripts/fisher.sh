#!/bin/bash
#PBS -N fisher
#PBS -l select=1:ncpus=50:ngpus=0:mem=64gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/1_coevolution
cd /home/lnemati/pathway_crosstalk/code/1_coevolution
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python fisher.py
exit 0
