#!/bin/bash
#PBS -N aggr
#PBS -l select=1:ncpus=2:ngpus=0:mem=32gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/6_gpcrs
cd /home/lnemati/pathway_crosstalk/code/6_gpcrs
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python aggregate_gpcrs.py
exit 0
