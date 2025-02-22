#!/bin/bash
#PBS -N surv_paraganglia
#PBS -l select=1:ncpus=1:ngpus=0:mem=15gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/7_survival
cd /home/lnemati/pathway_crosstalk/code/7_survival
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate r_survival
Rscript survival.R /home/lnemati/pathway_crosstalk/data/survival_data/tissues/paraganglia.csv /home/lnemati/pathway_crosstalk/results/survival/tissues
exit 0
