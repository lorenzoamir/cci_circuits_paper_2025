#!/bin/bash
#PBS -N pathways_cooccurrences
#PBS -l select=1:ncpus=8:ngpus=0:mem=16gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/1_compare
cd /home/lnemati/pathway_crosstalk/code/1_compare
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python pathways_cooccurrences.py
exit 0
