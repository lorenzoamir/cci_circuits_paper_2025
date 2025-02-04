#!/bin/bash
#PBS -N aggr
#PBS -l select=1:ncpus=2:ngpus=0:mem=8gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/7_survival
cd /home/lnemati/pathway_crosstalk/code/7_survival
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python aggregate.py --inputdir /home/lnemati/pathway_crosstalk/results/survival/tissues --outputdir /home/lnemati/pathway_crosstalk/results/survival/aggregate
exit 0
