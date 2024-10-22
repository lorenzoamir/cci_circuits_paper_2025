#!/bin/bash
#PBS -N wgcna_pipeline
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source "/projects/bioinformatics/snsutils/snsutils.sh"

SUBSET=0

SEPARATE=0 # Cant't run SEPARATE with other steps because it wont detect the .h5ad files

DESEQ=0
WGCNA=0
NETWORK=0
HUBS=0
RANK_ADJ=0
COEVOLUTION=0 # Check if this shuld be kept now that you have the new coev pipeline
INTERACTIONS=0
SCORE_INTS=1
ENRICHMENT=0
STATS=0

SEPARATE_QUEUE="q02anacreon"
DESEQ_QUEUE="q02anacreon"
WGCNA_QUEUE="q02anacreon"
NETWORK_QUEUE="q02gaia"
HUBS_QUEUE="q02anacreon"
RANK_ADJ_QUEUE="q02anacreon"
COEVOLUTION_QUEUE="q02anacreon"
INTERACTIONS_QUEUE="q02anacreon"
ENRICHMENT_QUEUE="q02gaia"
STATS_QUEUE="q02anacreon"

SEP_NCPUS=16
DESEQ_NCPUS=8
WGCNA_NCPUS=4
NETWORK_NCPUS=2
HUBS_NCPUS=2
RANK_ADJ_NCPUS=2
COEVOLUTION_NCPUS=8
INTERACTIONS_NCPUS=4
ENRICHMENT_NCPUS=4
STATS_NCPUS=8

SEP_MEMORY=64gb
DESEQ_MEMORY=10gb
WGCNA_MEMORY=16gb
NETWORK_MEMORY=32gb
HUBS_MEMORY=20gb
RANK_ADJ_MEMORY=16gb # Failed with 12gb
COEVOLUTION_MEMORY=18gb
INTERACTIONS_MEMORY=8gb
ENRICHMENT_MEMORY=6gb
STATS_MEMORY=24gb # Failed with 16gb, succeeded with 24gb

lr_resource="/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv"

cd "/home/lnemati/pathway_crosstalk/code/0_wgcna"

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

# ----- SEPARATE -----
if [ $SEPARATE -eq 1 ]; then
    echo "SEPARATE"

    # Create job script
    separate_name=separate
    separate_script=/home/lnemati/pathway_crosstalk/code/0_wgcna/separate_tissues.sh

    separate_id=$(fsub \
        -p "$separate_script" \
        -n "$separate_name" \
        -nc "$SEP_NCPUS" \
        -m "$SEP_MEMORY" \
        -q "$SEPARATE_QUEUE" \
        -e "WGCNA" \
        -c "python separate_tissues.py --inputfile /projects/bioinformatics/DB/Xena/TCGA_GTEX/TCGA_GTEX.h5ad --outputdir /projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/")

    echo ""
fi

# ----- Analysis -----

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA

# find all .h5ad files in the data directory and its subdirectories
data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
mapfile -t files < <(find "$data_dir" -name "*.h5ad")

# If SUBSET is positive, only use the first SUBSET files
if [ $SUBSET -gt 0 ]; then
    files=("${files[@]:0:$SUBSET}")
fi

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
 
        # Wait for separate job to finish
        waiting_list=""
        [ $SEPARATE -eq 1 ] && waiting_list="$waiting_list:$separate_id"
        echo "Waiting list: $waiting_list"

        deseq_id=$(fsub \
            -p "$deseq_script" \
            -n "$deseq_name" \
            -nc "$DESEQ_NCPUS" \
            -m "$DESEQ_MEMORY" \
            -q "$DESEQ_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python deseq_norm.py --input ${file} --ncpus $DESEQ_NCPUS")

        echo ""
    fi

    # ----- WGCNA -----
    if [ $WGCNA -eq 1 ]; then
        echo "WGCNA"

        # Create job script
        wgcna_name=wgcna_"$job_name"
        wgcna_script=$(dirname "$file")/scripts/$wgcna_name.sh
        
        # Wait for deseq job to finish
        waiting_list=""
        [ $DESEQ -eq 1 ] && waiting_list="$waiting_list:$deseq_id"
        echo "Waiting list: $waiting_list"

        wgcna_id=$(fsub \
            -p "$wgcna_script" \
            -n "$wgcna_name" \
            -nc "$WGCNA_NCPUS" \
            -m "$WGCNA_MEMORY" \
            -q "$WGCNA_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python wgcna.py --input ${file}")

        echo ""
    fi

    wgcnafile=$(dirname "$file")/wgcna_"$job_name".p

    # ----- NETWORK -----
    if [ $NETWORK -eq 1 ]; then
        echo "NETWORK"

        # Create job script
        network_name=network_"$job_name"
        network_script=$(dirname "$file")/scripts/$network_name.sh
        
        # Wait for wgcna job to finish
        waiting_list=""
        [ $WGCNA -eq 1 ] && waiting_list="$waiting_list:$wgcna_id"
        echo "Waiting list: $waiting_list"

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

    # ----- HUBS -----
    if [ $HUBS -eq 1 ]; then
        echo "HUBS"

        # Create job script
        hubs_name=hubs_"$job_name"
        hubs_script=$(dirname "$file")/scripts/$hubs_name.sh

        # Wait for wgcna job to finish
        waiting_list=""
        [ $WGCNA -eq 1 ] && waiting_list="$waiting_list:$wgcna_id"
        echo "Waiting list: $waiting_list"

        hubs_id=$(fsub \
            -p "$hubs_script" \
            -n "$hubs_name" \
            -nc "$HUBS_NCPUS" \
            -m "$HUBS_MEMORY" \
            -q "$HUBS_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python hubs_connectivities.py --input ${wgcnafile}")

        echo ""
    fi

    # ----- RANK ADJACENCY -----
    if [ $RANK_ADJ -eq 1 ]; then
        echo "RANK ADJACENCY"

        # Create job script
        rank_adj_name=rank_adj_"$job_name"
        rank_adj_script=$(dirname "$file")/scripts/$rank_adj_name.sh

        # Wait for wgcna job to finish
        waiting_list=""
        [ $WGCNA -eq 1 ] && waiting_list="$waiting_list:$wgcna_id"
        echo "Waiting list: $waiting_list"

        rank_adj_id=$(fsub \
            -p "$rank_adj_script" \
            -n "$rank_adj_name" \
            -nc "$RANK_ADJ_NCPUS" \
            -m "$RANK_ADJ_MEMORY" \
            -q "$RANK_ADJ_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python rank_pairs.py --input ${wgcnafile}")

        echo ""
    fi

    # ----- INTERACTIONS -----
    if [ $INTERACTIONS -eq 1 ]; then
        echo "INTERACTIONS"

        # Create job script
        interactions_name=interactions_"$job_name"
        interactions_script=$(dirname "$file")/scripts/$interactions_name.sh
       
        # Wait for wgcna job to finish
        waiting_list=""
        [ $WGCNA -eq 1 ] && waiting_list="$waiting_list:$wgcna_id"
        echo "Waiting list: $waiting_list"

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

    # ----- SCORE INTERACTIONS -----
    if [ $SCORE_INTS -eq 1 ]; then
        echo "SCORE INTERACTIONS"

        # Create job script
        score_ints_name=score_ints_"$job_name"
        score_ints_script=$(dirname "$file")/scripts/$score_ints_name.sh

        # Wait for wgcna and interactions jobs to finish
        waiting_list=""
        [ $WGCNA -eq 1 ] && waiting_list="$waiting_list:$wgcna_id"
        echo "Waiting list: $waiting_list"

        score_ints_id=$(fsub \
            -p "$score_ints_script" \
            -n "$score_ints_name" \
            -nc "$INTERACTIONS_NCPUS" \
            -m "$INTERACTIONS_MEMORY" \
            -q "$INTERACTIONS_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python score_interactions.py --input ${wgcnafile}")

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
     
        # Wait for wgcna job to finish
        waiting_list=""
        [ $WGCNA -eq 1 ] && waiting_list="$waiting_list:$wgcna_id"
        [ $INTERACTIONS -eq 1 ] && waiting_list="$waiting_list:$interactions_id" 
        echo "Waiting list: $waiting_list"

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

        # Wait for wgcna and interactions jobs to finish
        waiting_list=""
        [ $WGCNA -eq 1 ] && waiting_list="$waiting_list:$wgcna_id"
        [ $INTERACTIONS -eq 1 ] && waiting_list="$waiting_list:$interactions_id"
        echo "Waiting list: $waiting_list"

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

