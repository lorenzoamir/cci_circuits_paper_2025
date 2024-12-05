#!/bin/bash
#PBS -N all_head_and_neck_region
#PBS -l select=1:ncpus=2:ngpus=0:mem=8gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/7_survival
cd /home/lnemati/pathway_crosstalk/code/7_survival
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python survival_all.py --tissuefile /home/lnemati/pathway_crosstalk/data/survival_data/head_and_neck_region.csv
exit 0
