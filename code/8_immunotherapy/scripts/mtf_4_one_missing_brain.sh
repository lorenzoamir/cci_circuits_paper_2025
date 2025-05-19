#!/bin/bash
#PBS -N mtf_4_one_missing_brain
#PBS -l select=1:ncpus=1:ngpus=1:mem=16gb
#PBS -q q07gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/8_immunotherapy
cd /home/lnemati/pathway_crosstalk/code/8_immunotherapy
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate tabpfn
python individual.py --motif 4_one_missing --data /home/lnemati/pathway_crosstalk/data/immunotherapy/tissues/brain.csv --outdir /home/lnemati/pathway_crosstalk/results/immunotherapy/individual_interactions_and_motifs/
exit 0
