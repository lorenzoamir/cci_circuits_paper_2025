#!/bin/bash
#PBS -N dns_head_and_neck_region_tumor
#PBS -l select=1:ncpus=8:ngpus=0:mem=10gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
cd /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python density.py --motifs /home/lnemati/pathway_crosstalk/results/crosstalk/motifs_per_tissue/tumor/head_and_neck_region/motifs.csv --network /home/lnemati/pathway_crosstalk/data/networks/tumor/head_and_neck_region.csv.gz
exit 0
