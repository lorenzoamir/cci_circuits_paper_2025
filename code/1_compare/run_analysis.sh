#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

COMPARE=0

STATS=1
STATS_QUEUE='q02anacreon'
STATS_NCPUS=8
STATS_MEMORY=32gb

SANKEY=0
SANKEY_QUEUE='q02gaia'
SANKEY_NCPUS=4
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
    tumor_interactions=$(find "$directory" -name "LR_interactions.csv" -type f)

    # read normal file location from corresponding_normal_wgcna.txt
    # if corresponding_normal_wgcna.txt doesn't exist, print a message and skip
    normal_wgcna=$(cat "$file")
    normal_interactions=$(find "$(dirname "$normal_wgcna")" -name "LR_interactions.csv" -type f)

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
        job_name="$(basename "$directory")_stats"
        # job script in the same directory as the tumor file
        job_script="$directory/scripts/stats.sh"
        echo "Creating job script: $job_script"

        touch "$job_script"

        echo "#!/bin/bash" > "$job_script"
        echo "#PBS -l select=1:ncpus=$STATS_NCPUS:mem=$STATS_MEMORY" >> "$job_script"
        echo "#PBS -q $STATS_QUEUE" >> "$job_script"
        echo "#PBS -N $job_name" >> "$job_script"
        echo "" >> "$job_script"
        
        echo "eval \"\$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)\"" >> "$job_script"
        echo "conda activate WGCNA" >> "$job_script"
        echo "cd /home/lnemati/pathway_crosstalk/code/1_compare/" >> "$job_script"
        echo "python /home/lnemati/pathway_crosstalk/code/1_compare/stats.py --tumor "$tumor_wgcna" --normal "$normal_wgcna"" >> "$job_script"

        echo "exit 0" >> "$job_script"
        
        # submit job
        echo "Submitting job for $job_script"
        qsub "$job_script"
        echo ""
    fi

    # ----- SANKEY -----
    if [ $SANKEY -eq 1 ]; then
        echo 'Sankey'

        # create job script for each tumor
        job_name="$(basename "$directory")_sankey"
        # job script in the same directory as the tumor file
        job_script="$directory/scripts/sankey.sh"
        echo "Creating job script: $job_script"

        touch "$job_script"

        echo "#!/bin/bash" > "$job_script"
        echo "#PBS -l select=1:ncpus=$SANKEY_NCPUS:mem=$SANKEY_MEMORY" >> "$job_script"
        echo "#PBS -q $SANKEY_QUEUE" >> "$job_script"
        echo "#PBS -N $job_name" >> "$job_script"
        echo "" >> "$job_script"
        
        echo "eval \"\$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)\"" >> "$job_script"
        echo "conda activate WGCNA" >> "$job_script"
        echo "cd /home/lnemati/pathway_crosstalk/code/1_compare/" >> "$job_script"
        echo "python /home/lnemati/pathway_crosstalk/code/1_compare/sankey.py --tumor "$tumor_interactions" --normal "$normal_interactions"" >> "$job_script"

        echo "exit 0" >> "$job_script"
        
        # submit job
        echo "Submitting job for $job_script"
        qsub "$job_script"
        echo ""
    fi
done
