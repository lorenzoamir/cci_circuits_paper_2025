#!/bin/bash
#PBS -N wgcna_pipeline
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source "/projects/bioinformatics/snsutils/snsutils.sh"

SUBSET=0

MODULARITY=1

MODULARITY_QUEUE="q02anacreon"

MODULARITY_NCPUS=2

MODULARITY_MEMORY=20gb

lr_resource="/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv"

cd "/home/lnemati/pathway_crosstalk/code/4_consensus"

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA

# find all wgcna objects
data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
mapfile -t files < <(find "$data_dir" -name "wgcna_*.p")

# If SUBSET is positive, only use the first SUBSET files
if [ $SUBSET -gt 0 ]; then
    files=("${files[@]:0:$SUBSET}")
fi

for file in "${files[@]}"; do
    echo "$(dirname "$file")"
    waiting_list=""
    # clean job name from extension .p
    job_name=$(basename "$file" | sed 's/.p//')
    # create scripts directory in the same directory as the input file
    mkdir -p $(dirname "$file")/scripts
   
    # ----- MODULARITY -----
    if [ $MODULARITY -eq 1 ]; then
        echo "MODULARITY"

        # Create job script
        modularity_name=modularity_"$job_name"
        modularity_script=$(dirname "$file")/scripts/$modularity_name.sh
        
        # Wait for wgcna job to finish
        waiting_list=""
        echo "Waiting list: $waiting_list"

        modularity_id=$(fsub \
            -p "$modularity_script" \
            -n "$modularity_name" \
            -nc "$MODULARITY_NCPUS" \
            -m "$MODULARITY_MEMORY" \
            -q "$MODULARITY_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python modularity.py --input $file")

        echo ""
    fi
done
exit 0

