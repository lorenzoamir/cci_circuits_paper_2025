#!/bin/bash
#PBS -N sankey_adrenal_gland
#PBS -l select=1:ncpus=1:ngpus=0:mem=8gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/2_analysis
cd /home/lnemati/pathway_crosstalk/code/2_analysis
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python sankey.py --tumor /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/adrenal_gland/tumor/adrenocortical_cancer/interactions/ccc.csv --normal /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/adrenal_gland/normal/adrenal_gland/interactions/ccc.csv
exit 0
