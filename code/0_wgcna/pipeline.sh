#!/bin/bash
#PBS -N wgcna_pipeline
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source "/projects/bioinformatics/snsutils/snsutils.sh"

SUBSET=0

SEPARATE=0 # Cant't run SEPARATE with other steps because it wont detect the .h5ad files

DESEQ=0
WGCNA=0
NETWORK=1
COEVOLUTION=1
INTERACTIONS=1
ENRICHMENT=1
STATS=1

SEPARATE_QUEUE="q02anacreon"
DESEQ_QUEUE="q02anacreon"
WGCNA_QUEUE="q02anacreon"
NETWORK_QUEUE="q02anacreon"
COEVOLUTION_QUEUE="q02anacreon"
INTERACTIONS_QUEUE="q02anacreon"
ENRICHMENT_QUEUE="q02gaia"
STATS_QUEUE="q02anacreon"

SEP_NCPUS=16
DESEQ_NCPUS=8
WGCNA_NCPUS=4
NETWORK_NCPUS=8
COEVOLUTION_NCPUS=8
INTERACTIONS_NCPUS=4
ENRICHMENT_NCPUS=4
STATS_NCPUS=8

SEP_MEMORY=64gb
DESEQ_MEMORY=10gb
WGCNA_MEMORY=16gb
NETWORK_MEMORY=12gb
COEVOLUTION_MEMORY=18gb # Failed with 16gb (only testis)
INTERACTIONS_MEMORY=7gb
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

    # ----- COEVOLUTION -----
    if [ $COEVOLUTION -eq 1 ]; then
        echo "COEVOLUTION"

        # Create job script
        coevolution_name=coevolution_"$job_name"
        coevolution_script=$(dirname "$file")/scripts/$coevolution_name.sh

        # Wait for network job to finish
        waiting_list=""
        [ $NETWORK -eq 1 ] && waiting_list="$waiting_list:$network_id"
        echo "Waiting list: $waiting_list"

        directory=$(dirname "$file")
        coevolution_matrix="/home/lnemati/resources/coevolution/jaccard_genes.csv.gz"

        coevolution_id=$(fsub \
            -p "$coevolution_script" \
            -n "$coevolution_name" \
            -nc "$COEVOLUTION_NCPUS" \
            -m "$COEVOLUTION_MEMORY" \
            -q "$COEVOLUTION_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python coevolution.py --directory ${directory} --input ${coevolution_matrix}") 

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

