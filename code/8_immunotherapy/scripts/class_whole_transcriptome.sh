#!/bin/bash
#PBS -N class_whole_transcriptome
#PBS -l select=1:ncpus=1:ngpus=1:mem=4gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/8_immunotherapy
cd /home/lnemati/pathway_crosstalk/code/8_immunotherapy
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate tabpfn
python classify.py --motif whole_transcriptome
exit 0
