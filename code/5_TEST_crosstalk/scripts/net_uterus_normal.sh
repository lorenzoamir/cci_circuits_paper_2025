#!/bin/bash
#PBS -N net_uterus_normal
#PBS -l select=1:ncpus=1:ngpus=0:mem=12gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
cd /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python aggregate_networks.py --inputdir /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/uterus/normal
exit 0
