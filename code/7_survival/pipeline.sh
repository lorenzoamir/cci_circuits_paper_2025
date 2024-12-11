#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

SURV=1 # run survial analysis for all ccc and crosstalk pairs

SURV_QUEUE='q02anacreon'

SURV_NCPUS=1

SURV_MEMORY=6gb

cd /home/lnemati/pathway_crosstalk/code/7_survival
script_dir="/home/lnemati/pathway_crosstalk/code/7_survival/scripts"
if [ ! -d "$script_dir" ]; then
    mkdir "$script_dir"
fi

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

data_dir=/home/lnemati/pathway_crosstalk/data/survival_data
waiting_list=""

# Find all .csv files in the data directory
tissue_files=$(find $data_dir -name "*.csv")

# Loop over all tissue files
for tissue_file in $tissue_files; do
    echo "Tissue file: $tissue_file"

    # Get the tissue name from the file name
    tissue=$(basename $tissue_file | sed 's/.csv//')

    # create job script for each tissue
    surv_name="surv_$tissue"
    surv_script="$script_dir/$surv_name.sh"

    surv_id=$(fsub \
        -p "$surv_script" \
        -n "$surv_name" \
        -nc "$SURV_NCPUS" \
        -m "$SURV_MEMORY" \
        -e "r_survival" \
        -q "$SURV_QUEUE" \
        -w "$waiting_list" \
	-c "Rscript survival.R $tissue_file"
    )
done

echo "Done: pipeline.sh"
