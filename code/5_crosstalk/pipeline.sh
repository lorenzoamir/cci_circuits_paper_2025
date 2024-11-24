#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

COEXP=0 # make co-expression of complexes network
MOTIFS=1 # find motifs in the network

COEXP_QUEUE='q02anacreon'
MOTIFS_QUEUE='q02anacreon'

COEXP_NCPUS=32
MOTIFS_NCPUS=50

COEXP_MEMORY=100gb
MOTIFS_MEMORY=90gb

cd /home/lnemati/pathway_crosstalk/code/5_crosstalk
script_dir="/home/lnemati/pathway_crosstalk/code/5_crosstalk/scripts"
if [ ! -d "$script_dir" ]; then
    mkdir "$script_dir"
fi

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
waiting_list=""

# FLOW takes all tissues at once
if [ $COEXP -eq 1 ]; then
    echo 'Coexp'

    # create job script for all tissues
    coexp_name="coexp"
    coexp_script="$script_dir/$coexp_name.sh"

    coexp_id=$(fsub \
        -p "$coexp_script" \
        -n "$coexp_name" \
        -nc "$COEXP_NCPUS" \
        -m "$COEXP_MEMORY" \
        -e "WGCNA" \
        -q "$COEXP_QUEUE" \
        -c "python coexp.py" )
fi

if [ $MOTIFS -eq 1 ]; then
    echo 'Motifs'

    # create job script for all tissues
    motifs_name="motifs"
    motifs_script="$script_dir/$motifs_name.sh"
    
    # if coexp is 1, add coexp to the waiting list 
    [[ $COEXP -eq 1 ]] && waiting_list="$coexp_id"

    motifs_id=$(fsub \
        -p "$motifs_script" \
        -n "$motifs_name" \
        -nc "$MOTIFS_NCPUS" \
        -m "$MOTIFS_MEMORY" \
        -e "WGCNA" \
        -q "$MOTIFS_QUEUE" \
        -w "$waiting_list" \
        -c "python motifs.py" )
fi

echo "Done: pipeline.sh"
