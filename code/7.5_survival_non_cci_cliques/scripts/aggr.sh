#!/bin/bash
#PBS -N aggr
#PBS -l select=1:ncpus=2:ngpus=0:mem=20gb
#PBS -q q02anacreon
mkdir -p /home/lnemati/pathway_crosstalk/code/7.5_survival_non_cci_cliques
cd /home/lnemati/pathway_crosstalk/code/7.5_survival_non_cci_cliques
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate wgcna_v2
python aggregate.py --inputdir /home/lnemati/pathway_crosstalk/results/survival/non_cci_cliques/tissues --outputdir /home/lnemati/pathway_crosstalk/results/survival/non_cci_cliques/aggregate
exit 0
