#!/bin/bash
#PBS -N enrichment_tumor_median
#PBS -l select=1:ncpus=2:ngpus=0:mem=8gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/4_consensus
cd /home/lnemati/pathway_crosstalk/code/4_consensus
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python enrichment.py --input /home/lnemati/pathway_crosstalk/results/consensus_modules/tumor/median/consensus_modules.csv
exit 0
