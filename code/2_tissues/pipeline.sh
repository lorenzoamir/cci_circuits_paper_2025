#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

COMPARE=0

COMPARE_QUEUE='q02anacreon'

COMPARE_NCPUS=4

COMPARE_MEMORY=8gb

cd /home/lnemati/pathway_crosstalk/code/1_analysis

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

# find all corresponding_normal_wgcna.txt files in the data directory
data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
mapfile -t files < <(find "$data_dir" -name "corresponding_normal_file.txt" -type f)

for file in "${files[@]}"; do
    echo "$file"
    waiting_list=""
    # Read file content and check if it contains only one line
    lines=$(wc -l < "$file")
    if [ "$lines" -gt 1 ]; then
        echo "Multiple normal files found for $file"
        echo "Skipping..."
        echo "-------------------------"
        continue
    fi

    directory=$(dirname "$file")
    echo "Tumor directory: $directory"
    # tumor file is the .p file in the same directory as the corresponding_normal_wgcna.txt, with the name starting with WGCNA_ 
    tumor_wgcna=$(find "$directory" -name "WGCNA_*.p" -type f)
    tumor_interactions=$(find "$directory" -name "interactions.csv" -type f)

    # read normal file location from corresponding_normal_wgcna.txt
    # if corresponding_normal_wgcna.txt doesn't exist, print a message and skip
    normal_wgcna=$(cat "$file")
    normal_interactions=$(find "$(dirname "$normal_wgcna")" -name "interactions.csv" -type f)

    if [ -z "$normal_wgcna" ]; then
        echo "Normal file not found for $file"
        echo "Skipping..."
        echo "-------------------------"
        continue
    fi
    
	# ----- STATS -----
    if [ $STATS -eq 1 ]; then
        echo 'Stats'

        # create job script for each tumor
        stats_name="$(basename "$directory")_stats"
        stats_script="$directory/scripts/stats.sh"

        stats_id=$(fsub \
            -p "$stats_script" \
            -n "$stats_name" \
            -nc "$STATS_NCPUS" \
            -m "$STATS_MEMORY" \
            -e "WGCNA" \
            -q "$STATS_QUEUE" \
            -w "$waiting_list" \
            -c "python stats.py --tumor $tumor_wgcna --normal $normal_wgcna")
    fi

    # ----- SANKEY -----
    if [ $SANKEY -eq 1 ]; then
        echo 'Sankey'

        # create job script for each tumor
        sankey_name="$(basename "$directory")_sankey"
        sankey_script="$directory/scripts/sankey.sh"

        sanky_id=$(fsub \
            -p "$sankey_script" \
            -n "$sankey_name" \
            -nc "$SANKEY_NCPUS" \
            -m "$SANKEY_MEMORY" \
            -e "WGCNA" \
            -q "$SANKEY_QUEUE" \
            -w "$waiting_list" \
            -c "python sankey.py --tumor $tumor_interactions --normal $normal_interactions")
    fi
done
