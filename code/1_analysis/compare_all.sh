#!/bin/bash
#PBS -N compare_all
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -q q02anacreon

# Annotate all interactions

ncpus=20
memory=40gb
queue=q02anacreon

extract_field() {
    field="$1"
    filename="$2"
    grep "$field" "$filename" | cut -d ':' -f 2 | sed 's/^[[:space:]]*//'
}

cd "/home/lnemati/pathway_crosstalk/code/1_analysis"

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA

# find all corresponding_normal_file.txt files in the data directory
data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
mapfile -t files < <(find "$data_dir" -name "corresponding_normal_file.txt" -type f)
	
results_file="/home/lnemati/pathway_crosstalk/results/comparison_results.csv"

# if results_file exists, copy it to results_file.bak and delete results_file
if [ -f "$results_file" ]; then
	echo "Backing up $results_file to $results_file.bak"
	cp "$results_file" "$results_file.bak"
	rm "$results_file"
fi

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
	
	# tumor file is the .p file in the same directory as the corresponding_normal_file.txt, with the name starting with WGCNA_
	tumor_file=$(find "$directory" -name "WGCNA_*.p" -type f)

	# read normal file location from corresponding_normal_file.txt
	# if corresponding_normal_file.txt doesn't exist, print a message and skip
	normal_file=$(cat "$file")
	if [ -z "$normal_file" ]; then
		echo "Normal file not found for $file"
		echo "Skipping..."
		echo "-------------------------"
		continue
	fi
	
	echo "Reading tumor general_info file"
	# read general_info.txt file in the same directory as the tumor file
	tumor_general_info="$directory/general_info.txt"
	echo "Tumor General Info File: $tumor_general_info"
	tumor_n_modules=$(extract_field "n_modules" "$tumor_general_info")
	echo "n_modules: $tumor_n_modules"
	tumor_n_genes=$(extract_field "n_genes" "$tumor_general_info")
	echo "n_genes: $tumor_n_genes"
	tumor_total_gene_pairs=$(extract_field "total_gene_pairs" "$tumor_general_info")
	echo "total_gene_pairs: $tumor_total_gene_pairs"
	tumor_gene_pairs_same_module=$(extract_field "gene_pairs_same_module" "$tumor_general_info")
	echo "gene_pairs_same_module: $tumor_gene_pairs_same_module"
	tumor_total_LR_pairs=$(extract_field "total_LR_pairs" "$tumor_general_info")
	echo "total_LR_pairs: $tumor_total_LR_pairs"
	tumor_LR_pairs_same_module=$(extract_field "LR_pairs_same_module" "$tumor_general_info")
	echo "LR_pairs_same_module: $tumor_LR_pairs_same_module"
	tumor_LR_pairs_Z_score=$(extract_field "LR_pairs_Z_score" "$tumor_general_info")
	echo "LR_pairs_Z_score: $tumor_LR_pairs_Z_score"
	echo ""

	echo "Reading normal general_info file"
	# read general_info.txt file in the same directory as the normal file
	normal_general_info="$(dirname "$normal_file")/general_info.txt"
	echo "Normal General Info File: $normal_general_info"
	normal_n_modules=$(extract_field "n_modules" "$normal_general_info")
	echo "n_modules: $normal_n_modules"
	normal_n_genes=$(extract_field "n_genes" "$normal_general_info")
	echo "n_genes: $normal_n_genes"
	normal_total_gene_pairs=$(extract_field "total_gene_pairs" "$normal_general_info")
	echo "total_gene_pairs: $normal_total_gene_pairs"
	normal_gene_pairs_same_module=$(extract_field "gene_pairs_same_module" "$normal_general_info")
	echo "gene_pairs_same_module: $normal_gene_pairs_same_module"
	normal_total_LR_pairs=$(extract_field "total_LR_pairs" "$normal_general_info")
	echo "total_LR_pairs: $normal_total_LR_pairs"
	normal_LR_pairs_same_module=$(extract_field "LR_pairs_same_module" "$normal_general_info")
	echo "LR_pairs_same_module: $normal_LR_pairs_same_module"
	normal_LR_pairs_Z_score=$(extract_field "LR_pairs_Z_score" "$normal_general_info")
	echo "LR_pairs_Z_score: $normal_LR_pairs_Z_score"
	echo ""	

	# write the results to the file
	tumor_name=$(basename "$directory")
	
	# write header
	if [ ! -s "$results_file" ]; then
		echo "Name,Tumor_n_modules,Normal_n_modules,Tumor_n_genes,Normal_n_genes,Tumor_total_gene_pairs,Normal_total_gene_pairs,Tumor_gene_pairs_same_module,Normal_gene_pairs_same_module,Tumor_total_LR_pairs,Normal_total_LR_pairs,Tumor_LR_pairs_same_module,Normal_LR_pairs_same_module,Tumor_LR_pairs_Z_score,Normal_LR_pairs_Z_score" > "$results_file"
	fi 

	echo "$tumor_name,$tumor_n_modules,$normal_n_modules,$tumor_n_genes,$normal_n_genes,$tumor_total_gene_pairs,$normal_total_gene_pairs,$tumor_gene_pairs_same_module,$normal_gene_pairs_same_module,$tumor_total_LR_pairs,$normal_total_LR_pairs,$tumor_LR_pairs_same_module,$normal_LR_pairs_same_module,$tumor_LR_pairs_Z_score,$normal_LR_pairs_Z_score" >> "$results_file"

	echo "-------------------------"

done

## find all tom.csv.gz files in the data directory
#
#mapfile -t files < <(find "$data_dir" -name "tom.csv.gz" -type f)
#
## if path to file contains /normal/ add "GTEx_" to the beginning of the file name
## if path to file contains /tumor/ add "TCGA_" to the beginning of the file name
#
#TOM_dir="/projects/bioinformatics/DB/Xena/TCGA_GTEX/TOMs_GTEX"
#
## create TOM directory if it doesn't existt
#if [ ! -d "$TOM_dir" ]; then
#	mkdir "$TOM_dir"
#	# Create GTEx and TCGA directories
#fi
#
## total number of files
#total_files=${#files[@]}
#i=0
#for file in "${files[@]}"; do
#	echo "Processing $file"
#	# print the number of files processed so far and the total number of files
#	i=$((i+1))
#	echo "File $i of $total_files"
#	# get the directory of the file
#	directory=$(dirname "$file")
#	# get the name of the file
#	file_name=$(basename "$file")
#	# if the directory ends with /normal it's a "GTEx" file
#	if [[ "$directory" == *"/normal" ]]; then
#		# add the name of the second to last directory in the path to the file name
#		# ./stomach/normal/tom.csv.gz becomes stomach_tom.csv.gz
#		prefix=$(basename $(dirname "$directory"))
#		new_file_name="$prefix"_"$file_name"
#		# copy the file to the TOM directory
#		echo "Copying $file to $TOM_dir/$new_file_name"
#		cp "$file" "$TOM_dir/$new_file_name"
#		echo "-------------------------"
#	else
#		echo "Not a GTEx file"
#		echo "Skipping..."
#		echo "-------------------------"
#	
#	fi
#
#	echo ""
#
#done

echo "Done"

exit 0

