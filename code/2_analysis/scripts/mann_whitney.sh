#!/bin/bash
#PBS -N mann_whitney
#PBS -l select=1:ncpus=16:ngpus=0:mem=24gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/2_analysis
cd /home/lnemati/pathway_crosstalk/code/2_analysis
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python rank_interactions.py
exit 0
