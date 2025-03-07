#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

SPLIT=1 # Train / Test split
CLASS=1 # classify immunotherapy response
AGGR=1 # aggregate all the results

SPLIT_QUEUE='q02anacreon'
CLASS_QUEUE='q02gaia'
AGGR_QUEUE='q02anacreon'

SPLIT_NCPUS=1
CLASS_NCPUS=1
AGGR_NCPUS=2

SPLIT_MEMORY=24gb
CLASS_MEMORY=4gb
AGGR_MEMORY=20gb

cd /home/lnemati/pathway_crosstalk/code/8_immunotherapy
script_dir="/home/lnemati/pathway_crosstalk/code/8_immunotherapy/scripts"
result_dir="/home/lnemati/pathway_crosstalk/results/immunotherapy"
if [ ! -d "$result_dir" ]; then
    mkdir "$result_dir"
fi

if [ ! -d "$script_dir" ]; then
    mkdir "$script_dir"
fi

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

waiting_list=""

# Split the data
if [ $SPLIT -eq 1 ]; then
    split_script="$script_dir/split.sh"
    split_id=$(fsub \
        -p "$split_script" \
        -n "split" \
        -nc "$SPLIT_NCPUS" \
        -m "$SPLIT_MEMORY" \
        -e "tabpfn" \
        -q "$SPLIT_QUEUE" \
        -c "python split.py"
    )
    waiting_list=$split_id
fi

#motifs_list=(whole_transcriptome all_ccis 4_triangle_extra 4_path 4_no_crosstalk 4_one_missing 4_clique 4_cycle 3_clique 3_path individual_ccis)
motifs_list=(whole_transcriptome all_ccis individual_ccis 4_triangle_extra 4_path 4_no_crosstalk 4_one_missing 4_clique 4_cycle 3_clique 3_path)

class_ids=""
# Loop over all motifs
if [ $CLASS -eq 1 ]; then
    for motif in "${motifs_list[@]}"; do
        echo "Motif: $motif"
        # create job script for each motif
        class_name="class_$motif"
        class_script="$script_dir/$class_name.sh"

        class_id=$(fsub \
            -p "$class_script" \
            -n "$class_name" \
            -nc "$CLASS_NCPUS" \
            -ng 1 \
            -e "tabpfn" \
            -m "$CLASS_MEMORY" \
            -q "$CLASS_QUEUE" \
            -w "$waiting_list" \
            -c "python classify.py --motif $motif"
            )

        class_ids="$class_ids:$class_id"
    done
fi

# Add class_ids to waiting_list
waiting_list=$waiting_list$class_ids

if [ $AGGR -eq 1 ]; then
    aggregate_dir="$result_dir/aggregate"

    aggr_script="$script_dir/aggr.sh"
    aggr_id=$(fsub \
        -p "$aggr_script" \
        -n "aggr" \
        -nc "$AGGR_NCPUS" \
        -m "$AGGR_MEMORY" \
        -e "WGCNA" \
        -q "$AGGR_QUEUE" \
        -w "$waiting_list" \
        -c "python aggregate.py"
    )
    waiting_list=$waiting_list:$aggr_id
fi

echo "Done: pipeline.sh"

