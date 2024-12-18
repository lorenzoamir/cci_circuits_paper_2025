#!/bin/bash
#PBS -N surv_pan_cancer
#PBS -l select=1:ncpus=1:ngpus=0:mem=6gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/7_survival
cd /home/lnemati/pathway_crosstalk/code/7_survival
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate r_survival
Rscript survival.R /home/lnemati/pathway_crosstalk/data/survival_data/pan_cancer.csv
exit 0
