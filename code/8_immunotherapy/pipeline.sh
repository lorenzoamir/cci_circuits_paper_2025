#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

SPLIT=0 # Separate the different datasets
FEATURES=0 # classify immunotherapy response with full feature sets e.g. all CCIs
INDIVIDUAL=0 # classify immunotherapy response with individual interactions/motifs
AGGR=1 # aggregate all the results

SPLIT_QUEUE='q02anacreon'
FEATURES_QUEUE='q02gaia'
INDIVIDUAL_QUEUE='q02gaia'
AGGR_QUEUE='q02anacreon'

SPLIT_NCPUS=1
FEATURES_NCPUS=1
INDIVIDUAL_NCPUS=1
AGGR_NCPUS=2

SPLIT_MEMORY=24gb
FEATURES_MEMORY=8gb
INDIVIDUAL_MEMORY=8gb
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
        -e "WGCNA" \
        -q "$SPLIT_QUEUE" \
        -c "python split.py"
    )
    waiting_list=$split_id
fi

#features_list=(whole_transcriptome all_signatures all_ccis all_motifs all_cliques signatures)
features_list=(signatures)

# Cohorts
#cohorts_list=($(find /home/lnemati/pathway_crosstalk/data/immunotherapy/cohorts -name "*.csv"))
#outdir='/home/lnemati/pathway_crosstalk/results/immunotherapy/cohorts'

# Tissues
cohorts_list=($(find /home/lnemati/pathway_crosstalk/data/immunotherapy/tissues -name "*.csv"))

outdir='/home/lnemati/pathway_crosstalk/results/immunotherapy/tissues'

# Subset to only the first cohorts
#cohorts_list=("${cohorts_list[@]:0:3}")

feat_ids=""
# Loop over all features and cohorts
if [ $FEATURES -eq 1 ]; then
    for data in "${cohorts_list[@]}"; do
        for features in "${features_list[@]}"; do
            # create job script for each feature and cohort
            feat_name="ftr_${features}_$(basename $data .csv)"
            feat_script="$script_dir/$feat_name.sh"

            feat_id=$(fsub \
                -p "$feat_script" \
                -n "$feat_name" \
                -nc "$FEATURES_NCPUS" \
                -ng 1 \
                -e "tabpfn" \
                -m "$FEATURES_MEMORY" \
                -q "$FEATURES_QUEUE" \
                -w "$waiting_list" \
                -c "python features.py --features $features --data $data --outdir $outdir"
                )

            feat_ids=$feat_ids:$feat_id 
        done
    done

    # Add feat_ids to waiting_list
    waiting_list=$waiting_list$feat_ids

fi

#motifs_list=(4_clique 3_clique individual_ccis 4_no_crosstalk 4_triangle_extra 4_path 4_one_missing 4_cycle 3_path)
motifs_list=(random_pairs)

outdir='/home/lnemati/pathway_crosstalk/results/immunotherapy/individual_interactions_and_motifs/'

# Loop over all motifs and cohorts
mtf_ids=""
if [ $INDIVIDUAL -eq 1 ]; then
    for motifs in "${motifs_list[@]}"; do
        for data in "${cohorts_list[@]}"; do
            # create job script for each motif and cohort
            motif_name="mtf_${motifs}_$(basename $data .csv)"
            motif_script="$script_dir/$motif_name.sh"

            motif_id=$(fsub \
                -p "$motif_script" \
                -n "$motif_name" \
                -nc "$INDIVIDUAL_NCPUS" \
                -ng 1 \
                -e "tabpfn" \
                -m "$INDIVIDUAL_MEMORY" \
                -q "$INDIVIDUAL_QUEUE" \
                -w "$waiting_list" \
                -c "python individual.py --motif $motifs --data $data --outdir $outdir"
                )

            mtf_ids=$mtf_ids:$motif_id
        done
    done

    # Add mtf_ids to waiting_list
    waiting_list=$waiting_list$mtf_ids

fi

# Find directories with results in /home/lnemati/pathway_crosstalk/results/immunotherapy/tissues
tissues_list=$(find /home/lnemati/pathway_crosstalk/results/immunotherapy/individual_interactions_and_motifs -maxdepth 1 -mindepth 1 -type d)
echo "Tissues list: $tissues_list"

if [ $AGGR -eq 1 ]; then
    for dir in $tissues_list; do
        tissue_name=$(basename $dir)
        echo "Aggregating results for $tissue_name"

        aggregate_dir="$result_dir/aggregate"
        aggr_name="aggr_${tissue_name}"
        aggr_script="$script_dir/$aggr_name.sh"

        aggr_id=$(fsub \
            -p "$aggr_script" \
            -n "$aggr_name" \
            -nc "$AGGR_NCPUS" \
            -m "$AGGR_MEMORY" \
            -e "WGCNA" \
            -q "$AGGR_QUEUE" \
            -w "$waiting_list" \
            -c "python aggregate.py --inputdir $dir" \
        )
        
    done
fi

echo "Done: pipeline.sh"

