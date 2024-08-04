#!/bin/bash
#PBS -N lr_pairs
#PBS -l select=1:ncpus=4:ngpus=0:mem=16gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/1_analysis
cd /home/lnemati/pathway_crosstalk/code/1_analysis
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python lr_pairs.py
exit 0
