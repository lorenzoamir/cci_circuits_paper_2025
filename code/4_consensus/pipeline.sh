#!/bin/bash
#PBS -N wgcna_pipeline
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source "/projects/bioinformatics/snsutils/snsutils.sh"

CLUSTER=0
ENRICHMENT=0
MODULARITY=1

CLUSTER_QUEUE="q02anacreon"
ENRICHMENT_QUEUE="q02gaia"
MODULARITY_QUEUE="q02gaia"

CLUSTER_NCPUS=4
ENRICHMENT_NCPUS=2
MODULARITY_NCPUS=2

CLUSTER_MEMORY=32gb
ENRICHMENT_MEMORY=8gb
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

conditions=("normal" "tumor")

# Iterate over tumor normal conditions
for condition in "${conditions[@]}"; do
    echo "Condition: $condition"

    if [ $CLUSTER -eq 1 ]; then
        echo "CLUSTER"

        # Create job script
        cluster_name=cluster_"$condition"
        cluster_script="./scripts/$cluster_name.sh"
        inputfile="/home/lnemati/pathway_crosstalk/results/consensus_modules/$condition/median_tom.csv"
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
    fi

    if [ $ENRICHMENT -eq 1 ]; then
        echo "ENRICHMENT"

        # Create job script
        enrichment_name=enrichment_"$condition"
        enrichment_script=./scripts/$enrichment_name.sh

        # Wait for cluster job to finish
        [ $CLUSTER -eq 1 ] && waiting_list="$cluster_id"
        
        inputfile="/home/lnemati/pathway_crosstalk/results/consensus_modules/$condition/consensus_modules.csv"

        enrichment_id=$(fsub \
            -p "$enrichment_script" \
            -n "$enrichment_name" \
            -nc "$ENRICHMENT_NCPUS" \
            -m "$ENRICHMENT_MEMORY" \
            -q "$ENRICHMENT_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python enrichment.py --input $inputfile" 
        )

        echo ""
    fi


    for file in "${files[@]}"; do
        echo "$(dirname "$file")"
        # clean job name from extension .p and prefix wgcna_
        job_name=$(basename "$file" | sed 's/wgcna_//g' | sed 's/.p//g')

        # ----- MODULARITY -----
        if [ $MODULARITY -eq 1 ]; then
            echo "MODULARITY"

            # Create job script
            modularity_name=md_"$job_name"_"$condition"
            modularity_script=./scripts/$modularity_name.sh
            
            # Wait for cluster job to finish
            [ $CLUSTER -eq 1 ] && waiting_list="$cluster_id"

            modularity_id=$(fsub \
                -p "$modularity_script" \
                -n "$modularity_name" \
                -nc "$MODULARITY_NCPUS" \
                -m "$MODULARITY_MEMORY" \
                -q "$MODULARITY_QUEUE" \
                -e "WGCNA" \
                -w "$waiting_list" \
                -c "python modularity.py --input $file --condition $condition")

            echo ""
        fi
    done
done

echo "Done: pipeline.sh"

exit 0


