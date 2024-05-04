#!/bin/bash
#PBS -N wgcna_master
#PBS -l select=1:ncpus=16:mem=24gb
#PBS -q q02anacreon

source "/projects/bioinformatics/snsutils/snsutils.sh"

SEPARATE=0

DESEQ=0
WGCNA=0
NETWORK=1
INTERACTIONS=0
ENRICHMENT=0
STATS=0

DESEQ_QUEUE="q02anacreon"
WGCNA_QUEUE="q02anacreon"
NETWORK_QUEUE="q02anacreon"
INTERACTIONS_QUEUE="q02anacreon"
ENRICHMENT_QUEUE="q02gaia"
STATS_QUEUE="q02anacreon"

DESEQ_NCPUS=8
WGCNA_NCPUS=4
NETWORK_NCPUS=4
INTERACTIONS_NCPUS=4
ENRICHMENT_NCPUS=4
STATS_NCPUS=8

DESEQ_MEMORY=10gb
WGCNA_MEMORY=16gb
NETWORK_MEMORY=12gb # Failed with 6gb
INTERACTIONS_MEMORY=6gb
ENRICHMENT_MEMORY=6gb
STATS_MEMORY=25gb # Failed with 16gb

lr_resource="/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv"

cd "/home/lnemati/pathway_crosstalk/code/0_wgcna"

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

# ----- SEPARATE -----
if [ $SEPARATE -eq 1 ]; then
    conda activate WGCNA
    # separate data according to tissue 
    echo "Separating data according to tissue"
    echo ""
    python /home/lnemati/pathway_crosstalk/code/0_wgcna/separate_tissues.py --inputfile "/projects/bioinformatics/DB/Xena/TCGA_GTEX/TCGA_GTEX.h5ad" --outputdir "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/"

    echo "Done"
    echo ""
fi

# ----- Analysis -----

# find all .h5ad files in the data directory and its subdirectories
data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
mapfile -t files < <(find "$data_dir" -name "*.h5ad")

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
 
for file in "${files[@]}"; do
    echo "$(dirname "$file")"
    waiting_list=""
    # clean job name
    job_name=$(basename "$file" | sed 's/.h5ad//g')
    # create scripts directory in the same directory as the input file
    mkdir -p $(dirname "$file")/scripts
   
    # ----- DESEQ2 ----- 
    if [ $DESEQ -eq 1 ]; then
        echo "DESEQ2"

        # Create job script
        deseq_name=deseq_"$job_name"
        deseq_script=$(dirname "$file")/scripts/$deseq_name.sh
        
        deseq_id=$(fsub \
            -p "$deseq_script" \
            -n "$deseq_name" \
            -nc "$DESEQ_NCPUS" \
            -m "$DESEQ_MEMORY" \
            -q "$DESEQ_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python deseq_norm.py --input ${file} --ncpus $DESEQ_NCPUS")

        waiting_list="$waiting_list:$deseq_id"
        echo "Waiting list: $waiting_list"
        echo ""
    fi

    # ----- WGCNA -----
    if [ $WGCNA -eq 1 ]; then
        echo "WGCNA"

        # Create job script
        wgcna_name=wgcna_"$job_name"
        wgcna_script=$(dirname "$file")/scripts/$wgcna_name.sh

        wgcna_id=$(fsub \
            -p "$wgcna_script" \
            -n "$wgcna_name" \
            -nc "$WGCNA_NCPUS" \
            -m "$WGCNA_MEMORY" \
            -q "$WGCNA_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python wgcna.py --input ${file}")

        waiting_list="$waiting_list:$wgcna_id"
        echo "Waiting list: $waiting_list"
        echo ""
    fi

    wgcnafile=$(dirname "$file")/wgcna_"$job_name".p

    # ----- NETWORK -----
    if [ $NETWORK -eq 1 ]; then
        echo "NETWORK"

        # Create job script
        network_name=network_"$job_name"
        network_script=$(dirname "$file")/scripts/$network_name.sh

        network_id=$(fsub \
            -p "$network_script" \
            -n "$network_name" \
            -nc "$NETWORK_NCPUS" \
            -m "$NETWORK_MEMORY" \
            -q "$NETWORK_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python network.py --input ${wgcnafile}")

        echo ""
    fi

    # ----- INTERACTIONS -----
    if [ $INTERACTIONS -eq 1 ]; then
        echo "INTERACTIONS"

        # Create job script
        interactions_name=interactions_"$job_name"
        interactions_script=$(dirname "$file")/scripts/$interactions_name.sh
        
        interactions_id=$(fsub \
            -p "$interactions_script" \
            -n "$interactions_name" \
            -nc "$INTERACTIONS_NCPUS" \
            -m "$INTERACTIONS_MEMORY" \
            -q "$INTERACTIONS_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python interactions.py --input ${wgcnafile} --lr_resource ${lr_resource}")

        echo ""
    fi 

    # ----- Enrichment analysis -----
    if [ $ENRICHMENT -eq 1 ]; then
        echo "ENRICHMENT"
        
        # input file is the wgcna output file
        echo "Input file for enrichment: $wgcnafile"

        # Create job script
        enrichment_name=enrichment_"$job_name"
        enrichment_script=$(dirname "$file")/scripts/$enrichment_name.sh
      
        enrichment_id=$(fsub \
            -p "$enrichment_script" \
            -n "$enrichment_name" \
            -nc "$ENRICHMENT_NCPUS" \
            -m "$ENRICHMENT_MEMORY" \
            -q "$ENRICHMENT_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python enrichment.py --input ${wgcnafile}")
        
        echo ""
    fi

    # ----- STATS -----
    if [ $STATS -eq 1 ]; then
        echo "STATS"

        # Create job script
        stats_name=stats_"$job_name"
        stats_script=$(dirname "$file")/scripts/$stats_name.sh

        [ $INTERACTIONS -eq 1 ] && waiting_list="$waiting_list:$interactions_id"

        stats_id=$(fsub \
            -p "$stats_script" \
            -n "$stats_name" \
            -nc "$STATS_NCPUS" \
            -m "$STATS_MEMORY" \
            -q "$STATS_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python stats.py --input ${wgcnafile}")

        echo ""
    fi
done

echo "Done: all jobs submitted"

exit 0

