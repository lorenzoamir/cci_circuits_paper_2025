#!/bin/bash
#PBS -N motifs_brain_substantia_nigra_normal
#PBS -l select=1:ncpus=8:ngpus=0:mem=10gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
cd /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python motifs.py --inputfile /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/brain/normal/brain_substantia_nigra/interactions/all_ccc_complex_pairs.csv
exit 0
