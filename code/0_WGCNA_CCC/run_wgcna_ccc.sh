#!/bin/bash
#PBS -N wgcna_master
#PBS -l select=1:ncpus=32:mem=64gb
#PBS -q q02anacreon

# Sections of the script
SEPARATE=0
WGCNA=1
CCC=0

# get ncpus from NCPUS environment variable
ncpus=20
memory=25gb
queue="q02anacreon"

lr_resource="/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/consensus.csv"

cd "/home/lnemati/pathway_crosstalk/code/0_WGCNA_CCC"

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

# ----- SEPARATE -----
if [ $SEPARATE -eq 1 ]; then
    conda activate WGCNA
    # separate data according to tissue 
    echo "Separating data according to tissue"
    echo ""
    python /home/lnemati/pathway_crosstalk/code/0_WGCNA_CCC/SeparateTissuesConditions.py --inputfile "/projects/bioinformatics/DB/Xena/TCGA_GTEX/TCGA_GTEX.h5ad" --outputdir "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/"

    echo "Done"
    echo ""
fi

# ----- WGCNA & CCC -----

# find all .h5ad files in the data directory and its subdirectories
data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
mapfile -t files < <(find "$data_dir" -name "*.h5ad")

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
 
for file in "${files[@]}"; do
    echo "$(dirname "$file")"
    # clean job name
    job_name=$(basename "$file" | sed 's/.h5ad//g')
    
    # ----- WGCNA -----

    if [ $WGCNA -eq 1 ]; then
        echo "WGCNA"

        # Create job script
        wgcna_name=wgcna_"$job_name"
        wgcna_script=$(dirname "$file")/$wgcna_name.sh
        echo "Creating job script for $wgcna_script"
        touch "$wgcna_script"

        echo "#!/bin/bash" > "$wgcna_script"
        echo "#PBS -l select=1:ncpus=$ncpus:mem=$memory" >> "$wgcna_script"
        echo "#PBS -q $queue" >> "$wgcna_script"
        echo "#PBS -N $wgcna_name" >> "$wgcna_script"
        echo "" >> "$wgcna_script"
        echo 'eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"' >> "$wgcna_script"
        echo "conda activate WGCNA" >> "$wgcna_script"
        echo "python /home/lnemati/pathway_crosstalk/code/0_WGCNA_CCC/wgcna.py --input "${file}" --ncpus ${ncpus} --lr_resource "${lr_resource}"" >> "${wgcna_script}" 
        echo "exit 0" >> "$wgcna_script"

        # submit job
        echo "Submitting job for $wgcna_script"
        wgcna_id=$(qsub "$wgcna_script")
        echo ""
    fi

    # ----- CCC -----
    echo "CCC"

    if [ $CCC -eq 1 ]; then

        # Create job script
        ccc_name=ccc_"$job_name"
        ccc_script=$(dirname "$file")/$ccc_name.sh
        echo "Creating job script for $ccc_script"
        touch "$ccc_script"

        echo "#!/bin/bash" > "$ccc_script"
        echo "#PBS -l select=1:ncpus=1:mem=$memory" >> "$ccc_script"
        echo "#PBS -q $queue" >> "$ccc_script"
        echo "#PBS -N $ccc_name" >> "$ccc_script"
        echo "" >> "$ccc_script"
        echo 'eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"' >> "$ccc_script"
        echo "conda activate c2c" >> "$ccc_script"
        echo "python /home/lnemati/pathway_crosstalk/code/0_WGCNA_CCC/ccc.py --input "${file}"" >> "${ccc_script}"
        echo "exit 0" >> "$ccc_script"

        # submit job
        echo "Submitting job for $ccc_script"
        if [ $WGNA -eq 1 ]; then
            qsub -W depend=afterok:$wgcna_id "$ccc_script"
        else
            qsub "$ccc_script"
        fi
        echo ""
    fi
done

echo "Done"

exit 0

