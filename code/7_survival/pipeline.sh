#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

SURV=1 # run survial analysis for all ccc and crosstalk pairs
AGGR=1 # aggregate all the results

SURV_QUEUE='q02anacreon'
AGGR_QUEUE='q02anacreon'

SURV_NCPUS=1
AGGR_NCPUS=2
PLOT_NCPUS=1

SURV_MEMORY=15gb
AGGR_MEMORY=20gb

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

#if [ $PLOT -eq 1 ]; then
#    plot_script="$script_dir/plot.sh"
#    plot_id=$(fsub \
#        -p "$plot_script" \
#        -n "plot" \
#        -nc "$PLOT_NCPUS" \
#        -m "$PLOT_MEMORY" \
#        -e "WGCNA" \
#        -q "$PLOT_QUEUE" \
#        -w "$waiting_list" \
#        -c "python plot.py"
#    )
#    waiting_list=$waiting_list:$plot_id
#fi

echo "Done: pipeline.sh"
