#!/bin/bash
#PBS -N ftr_whole_transcriptome_full_dataset
#PBS -l select=1:ncpus=1:ngpus=1:mem=32gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/8_immunotherapy
cd /home/lnemati/pathway_crosstalk/code/8_immunotherapy
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate tabpfn
python features.py --features whole_transcriptome --data /home/lnemati/pathway_crosstalk/data/immunotherapy/cohorts/full_dataset.csv --outdir /home/lnemati/pathway_crosstalk/results/immunotherapy/features_sets/cohorts
exit 0
