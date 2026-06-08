#!/bin/bash
#PBS -N ftr_all_ccis_ureter_imvigor210
#PBS -l select=1:ncpus=1:ngpus=1:mem=16gb
#PBS -q q07gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/8_immunotherapy
cd /home/lnemati/pathway_crosstalk/code/8_immunotherapy
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate tabpfn
python features.py --features all_ccis --data /home/lnemati/pathway_crosstalk/data/immunotherapy/cohorts/ureter_imvigor210.csv --outdir /home/lnemati/pathway_crosstalk/results/immunotherapy/features_sets/cohorts
exit 0
