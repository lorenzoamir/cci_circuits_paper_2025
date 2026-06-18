#!/bin/bash
#PBS -N mtf_3_non_cci_clique_brain
#PBS -l select=1:ncpus=2:ngpus=0:mem=16gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/8.5_immunotherapy_non_cci_cliques
cd /home/lnemati/pathway_crosstalk/code/8.5_immunotherapy_non_cci_cliques
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate tabpfn
python individual.py --motif 3_non_cci_clique --data /home/lnemati/pathway_crosstalk/data/immunotherapy/tissues/brain.csv --outdir /home/lnemati/pathway_crosstalk/results/immunotherapy/non_cci_cliques/individual_interactions_and_motifs/
exit 0
