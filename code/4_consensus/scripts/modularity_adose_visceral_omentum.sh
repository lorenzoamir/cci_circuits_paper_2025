#!/bin/bash
#PBS -N modularity_adose_visceral_omentum
#PBS -l select=1:ncpus=2:ngpus=0:mem=40gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/4_consensus
cd /home/lnemati/pathway_crosstalk/code/4_consensus
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python modularity.py --input /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/adipose_tissue/normal/adipose_visceral_omentum/wgcna_adipose_visceral_omentum.p --condition tumor --quantile median
exit 0
