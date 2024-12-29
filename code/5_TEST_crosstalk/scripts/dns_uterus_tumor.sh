#!/bin/bash
#PBS -N dns_uterus_tumor
#PBS -l select=1:ncpus=8:ngpus=0:mem=10gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
cd /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python density.py --motifs /home/lnemati/pathway_crosstalk/results/crosstalk/motifs_per_tissue/tumor/uterus/motifs.csv --network /home/lnemati/pathway_crosstalk/data/networks/tumor/uterus.csv.gz
exit 0
