#!/bin/bash
#PBS -N btp_head_and_neck_squamous_cell_carcinoma
#PBS -l select=1:ncpus=2:ngpus=0:mem=8gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/3_flow
cd /home/lnemati/pathway_crosstalk/code/3_flow
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python bootstrap_interactions.py --input /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/head_and_neck_region/tumor/head_and_neck_squamous_cell_carcinoma/wgcna_head_and_neck_squamous_cell_carcinoma.p
exit 0
