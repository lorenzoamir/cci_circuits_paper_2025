#!/bin/bash
#PBS -N cov_stomach
#PBS -l select=1:ncpus=1:ngpus=0:mem=15gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/7.5_survival_non_cci_cliques
cd /home/lnemati/pathway_crosstalk/code/7.5_survival_non_cci_cliques
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate r_survival
Rscript covariates_survival.R /home/lnemati/pathway_crosstalk/data/survival_data/tissues/stomach.csv
exit 0
