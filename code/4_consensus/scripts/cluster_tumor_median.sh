#!/bin/bash
#PBS -N cluster_tumor_median
#PBS -l select=1:ncpus=4:ngpus=0:mem=32gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/4_consensus
cd /home/lnemati/pathway_crosstalk/code/4_consensus
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python cluster.py --inputdir /home/lnemati/pathway_crosstalk/results/consensus_modules/tumor/median
exit 0
