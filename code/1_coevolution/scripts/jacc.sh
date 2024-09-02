#!/bin/bash
#PBS -N jacc
#PBS -l select=1:ncpus=32:ngpus=0:mem=32gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/1_coevolution
cd /home/lnemati/pathway_crosstalk/code/1_coevolution
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python jaccard.py --inputfile /home/lnemati/pathway_crosstalk/data/tumor_coev/all_cancers/Patient/mutations_df.csv
exit 0
