#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

ALL=1 # run survial analysis for all ccc and crosstalk pairs

ALL_QUEUE='q02anacreon'

ALL_NCPUS=2

ALL_MEMORY=8gb

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
    all_name="all_$tissue"
    all_script="$script_dir/$all_name.sh"

    all_id=$(fsub \
        -p "$all_script" \
        -n "$all_name" \
        -nc "$ALL_NCPUS" \
        -m "$ALL_MEMORY" \
        -e "WGCNA" \
        -q "$ALL_QUEUE" \
        -w "$waiting_list" \
        -c "python survival_all.py --tissuefile $tissue_file"
    )
done

echo "Done: pipeline.sh"
