#!/bin/bash
#PBS -N mtf_individual_ccis_stomach
#PBS -l select=1:ncpus=2:ngpus=1:mem=16gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/8.5_immunotherapy_non_cci_cliques
cd /home/lnemati/pathway_crosstalk/code/8.5_immunotherapy_non_cci_cliques
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate tabpfn
python individual.py --motif individual_ccis --data /home/lnemati/pathway_crosstalk/data/immunotherapy/tissues/stomach.csv --outdir /home/lnemati/pathway_crosstalk/results/immunotherapy/non_cci_cliques/individual_interactions_and_motifs/
exit 0
