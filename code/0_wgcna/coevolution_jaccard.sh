#!/bin/bash
#PBS -N coevolution_jaccard
#PBS -l select=1:ncpus=8:ngpus=0:mem=16gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/0_wgcna
cd /home/lnemati/pathway_crosstalk/code/0_wgcna
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python coevolution_jaccard.py
exit 0
