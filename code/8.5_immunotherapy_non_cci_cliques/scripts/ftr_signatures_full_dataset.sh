#!/bin/bash
#PBS -N ftr_signatures_full_dataset
#PBS -l select=1:ncpus=1:ngpus=1:mem=32gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/8_immunotherapy
cd /home/lnemati/pathway_crosstalk/code/8_immunotherapy
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate tabpfn
python features.py --features signatures --data /home/lnemati/pathway_crosstalk/data/immunotherapy/tissues/full_dataset.csv --outdir /home/lnemati/pathway_crosstalk/results/immunotherapy/features_sets/tissues
exit 0
