#!/bin/bash
#PBS -N rewiring
#PBS -l select=1:ncpus=4:ngpus=0:mem=4gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/1_compare
cd /home/lnemati/pathway_crosstalk/code/1_compare
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python rewiring.py
exit 0
