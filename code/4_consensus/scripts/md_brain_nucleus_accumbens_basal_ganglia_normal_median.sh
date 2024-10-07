#!/bin/bash
#PBS -N md_brain_nucleus_accumbens_basal_ganglia_normal_median
#PBS -l select=1:ncpus=2:ngpus=0:mem=40gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/4_consensus
cd /home/lnemati/pathway_crosstalk/code/4_consensus
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python modularity.py --input /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/brain/normal/brain_nucleus_accumbens_basal_ganglia/wgcna_brain_nucleus_accumbens_basal_ganglia.p --condition normal --quantile median
exit 0
