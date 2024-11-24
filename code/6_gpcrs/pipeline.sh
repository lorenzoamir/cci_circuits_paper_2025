#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

AGGREGATE=1 # aggregate the gpcrs.csv tables

AGGREGATE_QUEUE='q02anacreon'

AGGREGATE_NCPUS=2

AGGREGATE_MEMORY=32gb

cd /home/lnemati/pathway_crosstalk/code/6_gpcrs
script_dir="/home/lnemati/pathway_crosstalk/code/6_gpcrs/scripts"
if [ ! -d "$script_dir" ]; then
    mkdir "$script_dir"
fi

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

waiting_list=""

if [ $AGGREGATE -eq 1 ]; then
    echo 'Aggregate'

    # create job script for all tissues
    aggregate_name="aggr"
    aggregate_script="$script_dir/$aggregate_name.sh"

    aggregate_id=$(fsub \
        -p "$aggregate_script" \
        -n "$aggregate_name" \
        -nc "$AGGREGATE_NCPUS" \
        -m "$AGGREGATE_MEMORY" \
        -e "WGCNA" \
        -q "$AGGREGATE_QUEUE" \
        -c "python aggregate_gpcrs.py" )
fi

echo "Done: pipeline.sh"
