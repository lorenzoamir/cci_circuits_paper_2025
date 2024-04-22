#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

COMPARE=0
STATS=0
SANKEY=1

STATS_QUEUE='q02anacreon'
SANKEY_QUEUE='q02gaia'

STATS_NCPUS=8
SANKEY_NCPUS=4

STATS_MEMORY=32gb
SANKEY_MEMORY=4gb

cd /home/lnemati/pathway_crosstalk/code/1_compare/

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

# ----- Compare ----- 

if [ $COMPARE -eq 1 ]; then
    echo 'Compare'
    qsub compare_all.sh
fi

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

        #echo "Creating job script: $stats_script"

        #touch "$stats_script"
        #echo "#!/bin/bash" > "$stats_script"
        #echo "#PBS -l select=1:ncpus=$STATS_NCPUS:mem=$STATS_MEMORY" >> "$stats_script"
        #echo "#PBS -q $STATS_QUEUE" >> "$stats_script"
        #echo "#PBS -N $stas_name" >> "$stats_script"
        #echo "" >> "$stats_script"
        #
        #echo "eval \"\$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)\"" >> "$stats_script"
        #echo "conda activate WGCNA" >> "$stats_script"
        #echo "cd /home/lnemati/pathway_crosstalk/code/1_compare/" >> "$stats_script"
        #echo "python /home/lnemati/pathway_crosstalk/code/1_compare/stats.py --tumor "$tumor_wgcna" --normal "$normal_wgcna"" >> "$stats_script"

        #echo "exit 0" >> "$stats_script"
        #
        ## submit job
        #echo "Submitting job for $stats_script"
        #qsub "$stats_script"
        #echo ""
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

        #echo "Creating job script: $sankey_script"

        #touch "$sankey_script"

        #echo "#!/bin/bash" > "$sankey_script"
        #echo "#PBS -l select=1:ncpus=$SANKEY_NCPUS:mem=$SANKEY_MEMORY" >> "$sankey_script"
        #echo "#PBS -q $SANKEY_QUEUE" >> "$sankey_script"
        #echo "#PBS -N $sankey_name" >> "$sankey_script"
        #echo "" >> "$sankey_script"
        #
        #echo "eval \"\$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)\"" >> "$sankey_script"
        #echo "conda activate WGCNA" >> "$sankey_script"
        #echo "cd /home/lnemati/pathway_crosstalk/code/1_compare/" >> "$sankey_script"
        #echo "python /home/lnemati/pathway_crosstalk/code/1_compare/sankey.py --tumor "$tumor_interactions" --normal "$normal_interactions"" >> "$sankey_script"

        #echo "exit 0" >> "$sankey_script"
        #
        ## submit job
        #echo "Submitting job for $sankey_script"
        #qsub "$sankey_script"
        #echo ""
    fi
done
