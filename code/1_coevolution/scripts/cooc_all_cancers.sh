#!/bin/bash
#PBS -N cooc_all_cancers
#PBS -l select=1:ncpus=8:ngpus=0:mem=16gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/1_coevolution
cd /home/lnemati/pathway_crosstalk/code/1_coevolution
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python cooccurrence.py --inputfile /home/lnemati/pathway_crosstalk/data/tumor_coev/all_cancers/Patient/mutations_df.csv
exit 0
