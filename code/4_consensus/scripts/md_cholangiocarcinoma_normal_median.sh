#!/bin/bash
#PBS -N md_cholangiocarcinoma_normal_median
#PBS -l select=1:ncpus=2:ngpus=0:mem=40gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/4_consensus
cd /home/lnemati/pathway_crosstalk/code/4_consensus
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python modularity.py --input /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/bile_duct/tumor/cholangiocarcinoma/wgcna_cholangiocarcinoma.p --condition normal --quantile median
exit 0
