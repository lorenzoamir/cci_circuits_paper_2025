#!/bin/bash
#PBS -N separate
#PBS -l select=1:ncpus=8:ngpus=0:mem=8gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/1_coevolution
cd /home/lnemati/pathway_crosstalk/code/1_coevolution
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python separate.py
exit 0
