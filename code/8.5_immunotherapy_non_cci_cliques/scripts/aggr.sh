#!/bin/bash
#PBS -N aggr
#PBS -l select=1:ncpus=2:ngpus=0:mem=20gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/8_immunotherapy
cd /home/lnemati/pathway_crosstalk/code/8_immunotherapy
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python aggregate.py --inputdir /home/lnemati/pathway_crosstalk/data/immunotherapy/tissues/skin.csv
exit 0
