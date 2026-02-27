#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

SURV=0 # run survial analysis for all ccc and crosstalk pairs
AGGR=0 # aggregate all the results
COV=1

SURV_QUEUE='q02anacreon'
AGGR_QUEUE='q02anacreon'
COV_QUEUE='q02anacreon'

SURV_NCPUS=1
AGGR_NCPUS=2
COV_NCPUS=1

SURV_MEMORY=15gb
AGGR_MEMORY=20gb
COV_MEMORY=15gb

cd /home/lnemati/pathway_crosstalk/code/7_survival
script_dir="/home/lnemati/pathway_crosstalk/code/7_survival/scripts"
result_dir="/home/lnemati/pathway_crosstalk/results/survival"

if [ ! -d "$result_dir" ]; then
    mkdir "$result_dir"
fi

if [ ! -d "$script_dir" ]; then
    mkdir "$script_dir"
fi

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

data_dir=/home/lnemati/pathway_crosstalk/data/survival_data/tissues
tissues_dir="$result_dir/tissues" # output dir

waiting_list=""

# Find all .csv files in the data directory
tissue_files=$(find $data_dir -type f -name "*.csv")

waiting_list=""

# Loop over all tissue files
if [ $SURV -eq 1 ]; then
    for tissue_file in $tissue_files; do
        echo "Tissue file: $tissue_file"
        # Get the tissue name from the file name
        tissue=$(basename $tissue_file | sed 's/.csv//')

        # create job script for each tissue
        surv_name="surv_$tissue"
        surv_script="$script_dir/$surv_name.sh"

        # if tissuefile is pan_cancer.csv, then increase memory
        actual_memory="$SURV_MEMORY"
        if [ "$tissue" == "pan_cancer" ]; then
            actual_memory="16gb"
        fi

        surv_id=$(fsub \
            -p "$surv_script" \
            -n "$surv_name" \
            -nc "$SURV_NCPUS" \
            -e "r_survival" \
            -m "$SURV_MEMORY" \
            -q "$SURV_QUEUE" \
            -c "Rscript survival.R $tissue_file $tissues_dir"
            )

        waiting_list=$waiting_list:$surv_id
    done
fi

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
        -c "python aggregate.py --inputdir $tissues_dir --outputdir $aggregate_dir"
    )
    waiting_list=$waiting_list:$aggr_id
fi

if [ $COV -eq 1 ]; then
    for tissue_file in $tissue_files; do
        tissue=$(basename $tissue_file | sed 's/.csv//')
        cov_name="cov_$tissue"
        cov_script="$script_dir/$cov_name.sh"
        
        # wait for aggregate job 
        [ $AGGR -eq 1 ] && waiting_list=$waiting_list:$aggr_id 

        cov_id=$(fsub \
            -p "$cov_script" \
            -n "$cov_name" \
            -nc "$COV_NCPUS" \
            -m "$COV_MEMORY" \
            -e "r_survival" \
            -q "$COV_QUEUE" \
            -w "$waiting_list" \
            -c "Rscript covariates_survival.R $tissue_file"
        )

        waiting_list=$waiting_list
    done
fi

echo "Done: pipeline.sh"
