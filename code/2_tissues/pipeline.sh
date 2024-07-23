#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

SANKEY=1

SANKEY_QUEUE='q02anacreon'

SANKEY_NCPUS=8

SANKEY_MEMORY=8gb

cd /home/lnemati/pathway_crosstalk/code/2_tissues

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

# find all corresponding_normal_wgcna.txt files in the data directory
data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/

# Find all tissue directories
mapfile -t tissuedirs < <(find "$data_dir" -maxdepth 1 -mindepth 1 -type d)

for tissuedir in "${tissuedirs[@]}"; do
    echo "$tissuedir"
    
    # Inside each tissue directory, there are tumor and normal directories
    tumordir="$tissuedir/tumor"
    normaldir="$tissuedir/normal"

    # Check if both tumor and normal directories contain one and only one directory each
    tumorcount=$(find "$tumordir" -maxdepth 1 -mindepth 1 -type d | wc -l)
    normalcount=$(find "$normaldir" -maxdepth 1 -mindepth 1 -type d | wc -l)

    # Also check that none of the directories they contain are empty
    t_empty=$(find "$tumordir" -maxdepth 1 -mindepth 1 -type d -empty | wc -l)
    n_empty=$(find "$normaldir" -maxdepth 1 -mindepth 1 -type d -empty | wc -l)

    if [ "$tumorcount" -eq 1 ] && [ "$normalcount" -eq 1 ] && [ "$t_empty" -eq 0 ] && [ "$n_empty" -eq 0 ]; then
        tumordir=$(find "$tumordir" -maxdepth 1 -mindepth 1 -type d)
        normaldir=$(find "$normaldir" -maxdepth 1 -mindepth 1 -type d)

        # Get name of normal directory
        name=$(basename "$normaldir")

    else
        echo "Skipping"
        continue
    fi

    # Run pipeline for each tissue
    if [ $SANKEY -eq 1 ]; then
        echo 'Sankey'

        # create job script for each tumor
        sankey_name="sankey_$name"
        sankey_script="$tumordir/scripts/sankey.sh"

        sankey_id=$(fsub \
            -p "$sankey_script" \
            -n "$sankey_name" \
            -nc "$SANKEY_NCPUS" \
            -m "$SANKEY_MEMORY" \
            -e "WGCNA" \
            -q "$SANKEY_QUEUE" \
            -c "python sankey.py --tumordir $tumordir --normaldir $normaldir")
    fi

done







# OLD CODE

: '

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
'
