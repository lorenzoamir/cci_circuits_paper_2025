#!/bin/bash
#PBS -N rank_rewired
#PBS -l select=1:ncpus=8:ngpus=0:mem=16gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/2_tissues
cd /home/lnemati/pathway_crosstalk/code/2_tissues
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python rank_rewired.py
exit 0
