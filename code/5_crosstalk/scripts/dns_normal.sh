#!/bin/bash
#PBS -N dns_normal
#PBS -l select=1:ncpus=16:ngpus=0:mem=10gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/5_crosstalk
cd /home/lnemati/pathway_crosstalk/code/5_crosstalk
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python density.py --motifs /home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/motifs/normal/motifs.csv
exit 0
