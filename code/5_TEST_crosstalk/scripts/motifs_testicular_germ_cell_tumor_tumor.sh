#!/bin/bash
#PBS -N motifs_testicular_germ_cell_tumor_tumor
#PBS -l select=1:ncpus=8:ngpus=0:mem=10gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
cd /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
python motifs.py --inputfile /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/testis/tumor/testicular_germ_cell_tumor/interactions/all_ccc_complex_pairs.csv
exit 0
