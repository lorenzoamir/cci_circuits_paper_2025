#!/bin/bash
#PBS -N separate
#PBS -l select=1:ncpus=16:ngpus=0:mem=64gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/0_wgcna
cd /home/lnemati/pathway_crosstalk/code/0_wgcna
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python separate_tissues.py --inputfile /projects/bioinformatics/DB/Xena/TCGA_GTEX/TCGA_GTEX.h5ad --outputdir /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
exit 0
