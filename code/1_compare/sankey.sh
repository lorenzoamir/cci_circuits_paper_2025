#!/bin/bash
#PBS -N sankey_master 
#PBS -l select=1:ncpus=8:mem=16gb
#PBS -q q02anacreon

# Annotate all interactions

ncpus=4
memory=8gb
queue=q02anacreon

cd "/home/lnemati/pathway_crosstalk/code/1_compare/"

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA

# find all corresponding_normal_file.txt files in the data directory
data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
mapfile -t files < <(find "$data_dir" -name "corresponding_normal_file.txt" -type f)

for file in "${files[@]}"; do
	
	# Read file content and check if it contains more than one line, if it does skip
	lines=$(wc -l < "$file")
	if [ "$lines" -gt 1 ]; then
		echo "Multiple normal files found for $file"
		echo "Skipping..."
		echo "-------------------------"
		continue
	fi

	directory=$(dirname "$file")
	echo "Tumor directory: $directory"
	
	# tumor file is the LR_interactions.csv file in the same directory as the corresponding_normal_file.txt
	tumor_file=$(find "$directory" -name "LR_interactions.csv" -type f)

	# read normal file location from corresponding_normal_file.txt
	# if corresponding_normal_file.txt doesn't exist, print a message and skip
    # else substitute the normal wgcna file with the LR_interactions.csv file in the same directory
	normal_file=$(cat "$file")
    normal_dir=$(dirname "$normal_file")
	if [ -z "$normal_file" ]; then
		echo "Normal file not found for $file"
		echo "Skipping..."
		echo "-------------------------"
		continue
    else
	    normal_file=$(find "$normal_dir" -name "LR_interactions.csv" -type f)
    fi
	
	# create job script for each tumor
	job_name="$(basename "$directory")_sankey"
	# job script in the same directory as the tumor file
	job_script="$directory/sankey.sh"
	echo "Creating job script: $job_script"

	touch "$job_script"

	echo "#!/bin/bash" > "$job_script"
	echo "#PBS -l select=1:ncpus=$ncpus:mem=$memory" >> "$job_script"
	echo "#PBS -q $queue" >> "$job_script"
	echo "#PBS -N $job_name" >> "$job_script"
	echo "" >> "$job_script"
	
	echo "eval \"\$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)\"" >> "$job_script"
	echo "conda activate WGCNA" >> "$job_script"
	echo "cd /home/lnemati/pathway_crosstalk/code/1_compare/" >> "$job_script"

	echo "python /home/lnemati/pathway_crosstalk/code/1_compare/sankey.py --tumor "$tumor_file" --normal "$normal_file"" >> "$job_script"

	echo "exit 0" >> "$job_script"
	# submit the job script	
	echo "Submitting job for $job_name"
	qsub "$job_script"
	echo ""

done

echo "Done"

exit 0

