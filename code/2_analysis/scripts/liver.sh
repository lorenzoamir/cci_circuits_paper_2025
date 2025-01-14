#!/bin/bash
#PBS -N liver
#PBS -l select=1:ncpus=1:ngpus=0:mem=8gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/2_analysis
cd /home/lnemati/pathway_crosstalk/code/2_analysis
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate base
python sankey.py --tumor /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/liver/tumor/liver_hepatocellular_carcinoma/interactions/ccc.csv --normal /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/liver/normal/liver/interactions/ccc.csv
exit 0
