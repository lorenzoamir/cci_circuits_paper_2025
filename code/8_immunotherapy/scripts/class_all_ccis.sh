#!/bin/bash
#PBS -N class_all_ccis
#PBS -l select=1:ncpus=1:ngpus=1:mem=8gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/8_immunotherapy
cd /home/lnemati/pathway_crosstalk/code/8_immunotherapy
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate tabpfn
python classify.py --motif all_ccis --data /home/lnemati/pathway_crosstalk/data/immunotherapy/cohorts/kidney_braun_2020.csv
exit 0
