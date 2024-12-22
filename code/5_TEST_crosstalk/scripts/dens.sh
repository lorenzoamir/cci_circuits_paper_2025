#!/bin/bash
#PBS -N dens
#PBS -l select=1:ncpus=24:ngpus=0:mem=32gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/5_crosstalk
cd /home/lnemati/pathway_crosstalk/code/5_crosstalk
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python density.py
exit 0
