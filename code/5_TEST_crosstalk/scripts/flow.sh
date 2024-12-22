#!/bin/bash
#PBS -N flow
#PBS -l select=1:ncpus=8:ngpus=0:mem=8gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/3_flow
cd /home/lnemati/pathway_crosstalk/code/3_flow
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python flow.py --tissues ""/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/uterus" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/kidney" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/brain" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/adrenal_gland" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/testis" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/liver" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/prostate" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/lung" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/breast" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/ovary" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/skin" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/pancreas" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/colon" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/esophagus" "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/stomach" "
exit 0
