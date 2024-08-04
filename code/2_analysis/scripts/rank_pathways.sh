#!/bin/bash
#PBS -N rank_pws
#PBS -l select=1:ncpus=8:ngpus=0:mem=8gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/1_analysis
cd /home/lnemati/pathway_crosstalk/code/1_analysis
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python rank_pathways.py
exit 0
