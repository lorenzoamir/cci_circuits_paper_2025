#!/bin/bash
#PBS -N cls_whole_transcriptome
#PBS -l select=1:ncpus=1:ngpus=1:mem=8gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/8_immunotherapy
cd /home/lnemati/pathway_crosstalk/code/8_immunotherapy
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate tabpfn
python classify.py --motif whole_transcriptome --data /home/lnemati/pathway_crosstalk/data/immunotherapy/tissues/skin.csv --outdir /home/lnemati/pathway_crosstalk/results/immunotherapy/tissues
exit 0
