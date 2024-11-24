#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

COMPARE=0
COEVOLUTION=0
NODES=0
AVG=1
INT_NETWORK=0
INT_PAIRS=0
RANK_INT=0
LR_PAIRS=0
PW_NETWORK=0
RANK_PWS=0
CONSESUS=0
CORR_CONS=0

COMPARE_QUEUE='q02anacreon'
COEVOLUTION_QUEUE='q02anacreon'
NODES_QUEUE='q02gaia'
AVG_QUEUE='q02anacreon'
INT_NETWORK_QUEUE='q02anacreon'
INT_PAIRS_QUEUE='q02anacreon'
RANK_INT_QUEUE='q02gaia'
LR_PAIRS_QUEUE='q02anacreon'
PW_NETWORK_QUEUE='q02anacreon'
RANK_PWS_QUEUE='q02anacreon'
CONSESUS_QUEUE='q02gaia'
CORR_CONS_QUEUE='q02gaia'

COMPARE_NCPUS=8
COEVOLUTION_NCPUS=50
NODES_NCPUS=8
AVG_NCPUS=4
INT_NETWORK_NCPUS=8
INT_PAIRS_NCPUS=40
RANK_INT_NCPUS=32
LR_PAIRS_NCPUS=4
PW_NETWORK_NCPUS=8
RANK_PWS_NCPUS=8
CONSESUS_NCPUS=8
CORR_CONS_NCPUS=8

COMPARE_MEMORY=8gb
COEVOLUTION_MEMORY=64gb # Generating the full tensor requires 150gb
NODES_MEMORY=8gb
AVG_MEMORY=8gb
INT_NETWORK_MEMORY=16gb
INT_PAIRS_MEMORY=24gb
RANK_INT_MEMORY=24gb #Succeded with 24gb
LR_PAIRS_MEMORY=16gb
PW_NETWORK_MEMORY=16gb
RANK_PWS_MEMORY=8gb
CONSESUS_MEMORY=350gb # Succeded with 350gb
CORR_CONS_MEMORY=350gb

cd /home/lnemati/pathway_crosstalk/code/2_analysis
if [ ! -d "/home/lnemati/pathway_crosstalk/code/2_analysis/scripts" ]; then
    mkdir /home/lnemati/pathway_crosstalk/code/2_analysis/scripts
fi

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
# find all directories containing files starting with wgcna_ and ending with .p (case insensitive) in the data_dir
directories=$(find "$data_dir" -name "wgcna_*.p" -type f | xargs -n1 dirname | sort -u)
echo "Found $(echo "$directories" | wc -l) directories"

# make comma delimited string of directories
all_directories=$(echo "$directories" | tr '\n' ',' | sed 's/,$//')

waiting_list=""

if [ $COMPARE -eq 1 ]; then
    echo 'Compare'
    compare_name="compare"
    compare_script="/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/compare.sh"

    compare_id=$(fsub \
        -p "$compare_script" \
        -n "$compare_name" \
        -nc "$COMPARE_NCPUS" \
        -m "$COMPARE_MEMORY" \
        -e "WGCNA" \
        -q "$COMPARE_QUEUE" \
        -c "python compare.py --dir-list $all_directories")
fi

if [ $COEVOLUTION -eq 1 ]; then
    echo 'Coevolution'
    coevolution_name="coevolution"
    coevolution_script="/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/coevolution.sh"

    coevolution_id=$(fsub \
        -p "$coevolution_script" \
        -n "$coevolution_name" \
        -nc "$COEVOLUTION_NCPUS" \
        -m "$COEVOLUTION_MEMORY" \
        -e "WGCNA" \
        -q "$COEVOLUTION_QUEUE" \
        -c "python coevolution.py --use-existing --dir-list "$all_directories"")
fi

if [ $NODES -eq 1 ]; then
    echo 'Hubs'
    nodes_name="nodes"
    nodes_script="/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/nodes.sh"

    nodes_id=$(fsub \
        -p "$nodes_script" \
        -n "$nodes_name" \
        -nc "$NODES_NCPUS" \
        -m "$NODES_MEMORY" \
        -e "WGCNA" \
        -q "$NODES_QUEUE" \
        -c "python node_metrics.py --dir-list $all_directories")
fi

if [ $AVG -eq 1 ]; then
    echo 'Average'
    avg_name="avg"
    avg_script="/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/average.sh"

    avg_id=$(fsub \
        -p "$avg_script" \
        -n "$avg_name" \
        -nc "$AVG_NCPUS" \
        -m "$AVG_MEMORY" \
        -e "WGCNA" \
        -q "$AVG_QUEUE" \
        -c "python averages.py --dir-list $all_directories")
fi

if [ $INT_NETWORK -eq 1 ]; then
    echo 'Interactions Network'
    int_network_name="int_network"
    int_network_script="/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/interactions_network.sh"

    int_network_id=$(fsub \
        -p "$int_network_script" \
        -n "$int_network_name" \
        -nc "$INT_NETWORK_NCPUS" \
        -m "$INT_NETWORK_MEMORY" \
        -e "WGCNA" \
        -q "$INT_NETWORK_QUEUE" \
        -c "python interactions_network.py")
fi

if [ $INT_PAIRS -eq 1 ]; then
    echo 'Interactions Pairs'
    int_pairs_name="int_pairs"
    int_pairs_script="/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/pairs_of_interactions.sh"

    int_pairs_id=$(fsub \
        -p "$int_pairs_script" \
        -n "$int_pairs_name" \
        -nc "$INT_PAIRS_NCPUS" \
        -m "$INT_PAIRS_MEMORY" \
        -e "WGCNA" \
        -q "$INT_PAIRS_QUEUE" \
        -c "python pairs_of_interactions.py")
fi

if [ $RANK_INT -eq 1 ]; then
    echo 'Rank Interactions'
    rank_int_name="rank_int"
    rank_int_script="/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/rank_interactions.sh"

    rank_int_id=$(fsub \
        -p "$rank_int_script" \
        -n "$rank_int_name" \
        -nc "$RANK_INT_NCPUS" \
        -m "$RANK_INT_MEMORY" \
        -e "WGCNA" \
        -q "$RANK_INT_QUEUE" \
        -c "python rank_interactions.py")
fi

if [ $LR_PAIRS -eq 1 ]; then
    echo 'LR Pairs'
    lr_pairs_name="lr_pairs"
    lr_pairs_script="/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/lr_pairs.sh"

    lr_pairs_id=$(fsub \
        -p "$lr_pairs_script" \
        -n "$lr_pairs_name" \
        -nc "$LR_PAIRS_NCPUS" \
        -m "$LR_PAIRS_MEMORY" \
        -e "WGCNA" \
        -q "$LR_PAIRS_QUEUE" \
        -c "python lr_pairs.py")
fi

if [ $PW_NETWORK -eq 1 ]; then
    echo 'Pathways Network'
    pw_network_name="pw_network"
    pw_network_script="/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/pathways_network.sh"

    waiting_list=""
    [ $INT_NETWORK -eq 1 ] && waiting_list="$waiting_list:$int_network_id"

    pw_network_id=$(fsub \
        -p "$pw_network_script" \
        -n "$pw_network_name" \
        -nc "$PW_NETWORK_NCPUS" \
        -m "$PW_NETWORK_MEMORY" \
        -e "WGCNA" \
        -q "$PW_NETWORK_QUEUE" \
        -w "$waiting_list" \
        -c "python pathways_network.py")
fi

if [ $RANK_PWS -eq 1 ]; then
    echo 'Rank Pathways'
    rank_pws_name="rank_pws"
    rank_pws_script="/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/rank_pathways.sh"

    waiting_list=""
    [ $INT_NETWORK -eq 1 ] && waiting_list="$waiting_list:$int_network_id"
    [ $PW_NETWORK -eq 1 ] && waiting_list="$waiting_list:$pw_network_id"
    
    rank_pws_id=$(fsub \
        -p "$rank_pws_script" \
        -n "$rank_pws_name" \
        -nc "$RANK_PWS_NCPUS" \
        -m "$RANK_PWS_MEMORY" \
        -e "WGCNA" \
        -q "$RANK_PWS_QUEUE" \
        -w "$waiting_list" \
        -c "python rank_pathways.py")
fi

if [ $CONSESUS -eq 1 ]; then
    echo 'Consensus'
    consesus_name="consensus"
    consesus_script="/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/consensus.sh"
    
    # instead of passing the directories, we pass the actual wgcna files
    # follow the same format as the directories
    all_files=$(find "$data_dir" -name "wgcna_*.p" -type f | sort -u | tr '\n' ',' | sed 's/,$//')

    # run the job twice, one for tumor and one for normal
    consensus_tumor_id=$(fsub \
        -p "$consesus_script" \
        -n "$consesus_name"_tumor \
        -nc "$CONSESUS_NCPUS" \
        -m "$CONSESUS_MEMORY" \
        -e "WGCNA" \
        -q "$CONSESUS_QUEUE" \
        -c "python consensus_distance.py --condition tumor --file-list $all_files")

    consensus_normal_id=$(fsub \
        -p "$consesus_script" \
        -n "$consesus_name"_normal \
        -nc "$CONSESUS_NCPUS" \
        -m "$CONSESUS_MEMORY" \
        -e "WGCNA" \
        -q "$CONSESUS_QUEUE" \
        -c "python consensus_distance.py --condition normal --file-list $all_files")

fi

if [ $CORR_CONS -eq 1 ]; then
    echo 'Correlation Consensus'
    corr_cons_name="corr_cons"
    corr_cons_script="/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/${consesus_name}.sh"
    
    # instead of passing the directories, we pass the actual wgcna files
    # follow the same format as the directories
    all_files=$(find "$data_dir" -name "wgcna_*.p" -type f | sort -u | tr '\n' ',' | sed 's/,$//')

    # run the job twice, one for tumor and one for normal
    corr_cons_tumor_id=$(fsub \
        -p "$corr_cons_script" \
        -n "$corr_cons_name"_tumor \
        -nc "$CORR_CONS_NCPUS" \
        -m "$CORR_CONS_MEMORY" \
        -e "WGCNA" \
        -q "$CORR_CONS_QUEUE" \
        -c "python correlation_consensus.py --condition tumor --file-list $all_files")

    corr_cons_normal_id=$(fsub \
        -p "$corr_cons_script" \
        -n "$corr_cons_name"_normal \
        -nc "$CORR_CONS_NCPUS" \
        -m "$CORR_CONS_MEMORY" \
        -e "WGCNA" \
        -q "$CORR_CONS_QUEUE" \
        -c "python correlation_consensus.py --condition normal --file-list $all_files")
   
fi
echo "Done: all jobs submitted"

exit 0
