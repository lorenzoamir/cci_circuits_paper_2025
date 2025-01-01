#!/bin/bash
#PBS -N srv_brain_tumor
#PBS -l select=1:ncpus=1:ngpus=0:mem=6gb
#PBS -q q07anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
cd /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate r_survival
Rscript survival.R /home/lnemati/pathway_crosstalk/data/survival_data/brain.csv /home/lnemati/pathway_crosstalk/results/crosstalk/motifs_per_tissue/tumor/brain/motifs.csv
exit 0
