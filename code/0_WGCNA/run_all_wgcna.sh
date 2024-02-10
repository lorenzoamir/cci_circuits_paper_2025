#!/bin/bash
#PBS -N wgcna_master
#PBS -l select=1:ncpus=32:mem=64gb
#PBS -q q02anacreon

# Divide data according to tissue and run DeSeq2 normalization

# get ncpus from NCPUS environment variable
ncpus=20
memory=32gb
queue="q02anacreon"

cd "/home/lnemati/pathway_crosstalk/code/0_WGCNA"

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA

# separate data according to tissue 
echo "Separating data according to tissue"
echo ""
python /home/lnemati/pathway_crosstalk/code/0_WGCNA/SeparateTissuesConditions.py --inputfile "/projects/bioinformatics/DB/Xena/TCGA_GTEX/TCGA_GTEX.h5ad" --outputdir "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/"

echo "Done"
echo ""

# find all .h5ad files in the data directory and its subdirectories
data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
mapfile -t files < <(find "$data_dir" -name "*.h5ad")

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
 
for file in "${files[@]}"; do
	echo "$(dirname "$file")"
	# clean job name
	job_name=$(basename "$file" | sed 's/.h5ad//g')
    job_name=wgcna_"$job_name"
    # create a file with pbs job submission script for each file in the same folder as the file
	# job script will be named the same as the job with .sh extension (absolute path)
	job_script=$(dirname "$file")/$job_name.sh
	echo "Creating job script for $job_script"
	touch "$job_script"

	echo "#!/bin/bash" > "$job_script"
	echo "#PBS -l select=1:ncpus=$ncpus:mem=$memory" >> "$job_script"
	#echo "#PBS -l walltime=48:00:00" >> "$job_script"
	echo "#PBS -q $queue" >> "$job_script"
	echo "#PBS -N $job_name" >> "$job_script"
	echo "" >> "$job_script"
	echo 'eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"' >> "$job_script"
	echo "conda activate WGCNA" >> "$job_script"
	echo "python /home/lnemati/pathway_crosstalk/code/0_WGCNA/WGCNA.py --input "${file}" --ncpus ${ncpus}" >> "${job_script}" 
	echo "exit 0" >> "$job_script"
    # submit the job
	echo "Submitting job for $job_name"
    qsub "$job_script"
	echo ""
done

echo "Done"

exit 0

