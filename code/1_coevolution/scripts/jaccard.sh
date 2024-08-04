#!/bin/bash
#PBS -N jaccard
#PBS -l select=1:ncpus=8:ngpus=0:mem=18gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/1_coevolution
cd /home/lnemati/pathway_crosstalk/code/1_coevolution
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python coevolution_jaccard.py
exit 0
