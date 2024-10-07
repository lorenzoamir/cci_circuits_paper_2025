#!/bin/bash
#PBS -N md_skin_sun_eosed_lower_leg_normal_perc25
#PBS -l select=1:ncpus=2:ngpus=0:mem=40gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/4_consensus
cd /home/lnemati/pathway_crosstalk/code/4_consensus
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python modularity.py --input /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/skin/normal/skin_sun_exposed_lower_leg/wgcna_skin_sun_exposed_lower_leg.p --condition normal --quantile perc25
exit 0
