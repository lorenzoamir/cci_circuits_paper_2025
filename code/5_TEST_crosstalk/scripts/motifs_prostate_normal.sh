#!/bin/bash
#PBS -N motifs_prostate_normal
#PBS -l select=1:ncpus=8:ngpus=0:mem=10gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
cd /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python motifs.py --inputfile /home/lnemati/pathway_crosstalk/data/networks/normal/prostate.csv.gz
exit 0
