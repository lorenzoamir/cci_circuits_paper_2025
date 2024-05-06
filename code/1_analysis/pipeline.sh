#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

COMPARE=0
HUBS=0
INT_NETWORK=0
PW_NETWORK=0
RANK_PWS=1

COMPARE_QUEUE='q02anacreon'
HUBS_QUEUE='q02anacreon'
INT_NETWORK_QUEUE='q02anacreon'
PW_NETWORK_QUEUE='q02anacreon'
RANK_PWS_QUEUE='q02anacreon'

COMPARE_NCPUS=8
HUBS_NCPUS=8
INT_NETWORK_NCPUS=8
PW_NETWORK_NCPUS=8
RANK_PWS_NCPUS=8

COMPARE_MEMORY=8gb
HUBS_MEMORY=8gb
INT_NETWORK_MEMORY=16gb
PW_NETWORK_MEMORY=16gb
RANK_PWS_MEMORY=8gb

cd /home/lnemati/pathway_crosstalk/code/1_analysis
if [ ! -d "/home/lnemati/pathway_crosstalk/code/1_analysis/scripts" ]; then
    mkdir /home/lnemati/pathway_crosstalk/code/1_analysis/scripts
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
    # create job script for each tumor
    compare_name="compare"
    compare_script="/home/lnemati/pathway_crosstalk/code/1_analysis/scripts/compare.sh"

    compare_id=$(fsub \
        -p "$compare_script" \
        -n "$compare_name" \
        -nc "$COMPARE_NCPUS" \
        -m "$COMPARE_MEMORY" \
        -e "WGCNA" \
        -q "$COMPARE_QUEUE" \
        -c "python compare.py --dir_list $all_directories")
fi

if [ $HUBS -eq 1 ]; then
    echo 'Hubs'
    # create job script for each tumor
    hubs_name="hubs"
    hubs_script="/home/lnemati/pathway_crosstalk/code/1_analysis/scripts/hubs.sh"

    hubs_id=$(fsub \
        -p "$hubs_script" \
        -n "$hubs_name" \
        -nc "$HUBS_NCPUS" \
        -m "$HUBS_MEMORY" \
        -e "WGCNA" \
        -q "$HUBS_QUEUE" \
        -c "python hubs.py --dir_list $all_directories")
fi

if [ $INT_NETWORK -eq 1 ]; then
    echo 'Interactions Network'
    # create job script for each tumor
    int_network_name="int_network"
    int_network_script="/home/lnemati/pathway_crosstalk/code/1_analysis/scripts/interactions_network.sh"

    int_network_id=$(fsub \
        -p "$int_network_script" \
        -n "$int_network_name" \
        -nc "$INT_NETWORK_NCPUS" \
        -m "$INT_NETWORK_MEMORY" \
        -e "WGCNA" \
        -q "$INT_NETWORK_QUEUE" \
        -c "python interactions_network.py")
fi

if [ $PW_NETWORK -eq 1 ]; then
    echo 'Pathways Network'
    # create job script for each tumor
    pw_network_name="pw_network"
    pw_network_script="/home/lnemati/pathway_crosstalk/code/1_analysis/scripts/pathways_network.sh"

    # Wait for interactions network to finish
    waiting_list="$waiting_list:$int_network_id"

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
    # create job script for each tumor
    rank_pws_name="rank_pws"
    rank_pws_script="/home/lnemati/pathway_crosstalk/code/1_analysis/scripts/rank_pathways.sh"

    # Wait for pathways network to finish
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

exit 0


##mapfile -t files < <(find "$data_dir" -name "WGCNA_*.p" -type f)
#
#
#for file in "${files[@]}"; do
#    echo "$file"
#    waiting_list=""
#    # Read file content and check if it contains only one line
#    lines=$(wc -l < "$file")
#    if [ "$lines" -gt 1 ]; then
#        echo "Multiple normal files found for $file"
#        echo "Skipping..."
#        echo "-------------------------"
#        continue
#    fi
#
#    directory=$(dirname "$file")
#    echo "Tumor directory: $directory"
#    tumor_wgcna=$(find "$directory" -name "WGCNA_*.p" -type f)
#    tumor_interactions=$(find "$directory" -name "interactions.csv" -type f)
#
#    # if corresponding_normal_wgcna.txt doesn't exist, print a message and skip
#    normal_wgcna=$(cat "$file")
#    normal_interactions=$(find "$(dirname "$normal_wgcna")" -name "interactions.csv" -type f)
#
#    if [ -z "$normal_wgcna" ]; then
#        echo "Normal file not found for $file"
#        echo "Skipping..."
#        echo "-------------------------"
#        continue
#    fi
#    
#	# ----- STATS -----
#    if [ $STATS -eq 1 ]; then
#        echo 'Stats'
#
#        # create job script for each tumor
#        stats_name="$(basename "$directory")_stats"
#        stats_script="$directory/scripts/stats.sh"
#
#        stats_id=$(fsub \
#            -p "$stats_script" \
#            -n "$stats_name" \
#            -nc "$STATS_NCPUS" \
#            -m "$STATS_MEMORY" \
#            -e "WGCNA" \
#            -q "$STATS_QUEUE" \
#            -w "$waiting_list" \
#            -c "python stats.py --tumor $tumor_wgcna --normal $normal_wgcna")
#    fi
#
#    # ----- SANKEY -----
#    if [ $SANKEY -eq 1 ]; then
#        echo 'Sankey'
#
#        # create job script for each tumor
#        sankey_name="$(basename "$directory")_sankey"
#        sankey_script="$directory/scripts/sankey.sh"
#
#        sanky_id=$(fsub \
#            -p "$sankey_script" \
#            -n "$sankey_name" \
#            -nc "$SANKEY_NCPUS" \
#            -m "$SANKEY_MEMORY" \
#            -e "WGCNA" \
#            -q "$SANKEY_QUEUE" \
#            -w "$waiting_list" \
#            -c "python sankey.py --tumor $tumor_interactions --normal $normal_interactions")
#    fi
#done
