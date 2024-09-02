#!/bin/bash
#PBS -N dir_glioma
#PBS -l select=1:ncpus=8:ngpus=0:mem=16gb
#PBS -q q02gaia
mkdir -p /home/lnemati/pathway_crosstalk/code/1_coevolution
cd /home/lnemati/pathway_crosstalk/code/1_coevolution
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python directionality.py --mutationsdir /home/lnemati/pathway_crosstalk/data/tumor_coev/glioma --significant /home/lnemati/pathway_crosstalk/data/tumor_coev/glioma/Patient/cooccurrences_filtered.csv
exit 0
