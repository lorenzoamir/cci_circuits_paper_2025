#!/bin/bash
#PBS -N md_small_intestine_terminal_ileum_normal_perc25
#PBS -l select=1:ncpus=2:ngpus=0:mem=40gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/4_consensus
cd /home/lnemati/pathway_crosstalk/code/4_consensus
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python modularity.py --input /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/small_intestine/normal/small_intestine_terminal_ileum/wgcna_small_intestine_terminal_ileum.p --condition normal --quantile perc25
exit 0
