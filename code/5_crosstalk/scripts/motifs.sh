#!/bin/bash
#PBS -N motifs
#PBS -l select=1:ncpus=50:ngpus=0:mem=90gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/5_crosstalk
cd /home/lnemati/pathway_crosstalk/code/5_crosstalk
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python motifs.py
exit 0
