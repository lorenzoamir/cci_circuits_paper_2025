#!/bin/bash
#PBS -N wgcna_pipeline
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source "/projects/bioinformatics/snsutils/snsutils.sh"

SUBSET=0

CLUSTER=1
MODULARITY=0

CLUSTER_QUEUE="q02anacreon"
MODULARITY_QUEUE="q02gaia"

CLUSTER_NCPUS=4
MODULARITY_NCPUS=2

CLUSTER_MEMORY=32gb
MODULARITY_MEMORY=40gb

lr_resource="/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv"

cd "/home/lnemati/pathway_crosstalk/code/4_consensus"

# Make scripts directory if it does not exist
mkdir -p scripts

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA

# find all wgcna objects
data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
mapfile -t files < <(find "$data_dir" -name "wgcna_*.p")

# If SUBSET is positive, only use the first SUBSET files
if [ $SUBSET -gt 0 ]; then
    files=("${files[@]:0:$SUBSET}")
fi

conditions=("normal" "tumor")
quantiles=("perc25" "median")

if [ $CLUSTER -eq 1 ]; then
    echo "CLUSTER"

    # Iterate over the 4 possible combinations of conditions and quantiles
    for condition in "${conditions[@]}"; do
        for quantile in "${quantiles[@]}"; do
            echo "Condition: $condition"
            echo "Quantile: $quantile"

            # Create job script
            cluster_name=cluster_"$condition"_"$quantile"
            cluster_script="./scripts/$cluster_name.sh"
            inputdir="/home/lnemati/pathway_crosstalk/results/consensus_modules/$condition/$quantile"
            echo "$inputdir"

            cluster_id=$(fsub \
                -p "$cluster_script" \
                -n "$cluster_name" \
                -nc "$CLUSTER_NCPUS" \
                -m "$CLUSTER_MEMORY" \
                -q "$CLUSTER_QUEUE" \
                -e "WGCNA" \
                -w "$waiting_list" \
                -c "python cluster.py --inputdir $inputdir"
            )

            echo ""
        done
    done
fi

#  From now on we'll have to iterate over the two quantiles: median and perc25
quantiles=("perc25" "median")


for file in "${files[@]}"; do
    echo "$(dirname "$file")"
    waiting_list=""
    # clean job name from extension .p and prefix wgcna_
    job_name=$(basename "$file" | sed 's/wgcna_//g' | sed 's/.p//g')
    # create scripts directory in the same directory as the input file if it does not exist
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

echo "Done: pipeline.sh"

exit 0


