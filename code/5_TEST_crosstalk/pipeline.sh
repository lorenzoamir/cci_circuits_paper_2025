#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

NETWORKS=0 # aggregate subtissue networks into tissue networks
MOTIFS=1 # detect motifs in each tissue
#COEXP=0 # make co-expression of complexes network
#MOTIFS=0 # find motifs in the network
#DENS=0 # density

NETWORKS_QUEUE='q02anacreon'
MOTIFS_QUEUE='q02anacreon'
#COEXP_QUEUE='q02anacreon'
#MOTIFS_QUEUE='q02anacreon'
#DENS_QUEUE='q02anacreon'

NETWORKS_NCPUS=2
MOTIFS_NCPUS=8
#COEXP_NCPUS=32
#MOTIFS_NCPUS=50
#DENS_NCPUS=24

NETWORKS_MEMORY=16gb
MOTIFS_MEMORY=10gb
#COEXP_MEMORY=100gb
#MOTIFS_MEMORY=90gb
#DENS_MEMORY=32gb

cd /home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk/
script_dir="/home/lnemati/pathway_crosstalk/code/5_TEST_crosstalk/scripts"
if [ ! -d "$script_dir" ]; then
    mkdir "$script_dir"
fi

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
waiting_list=""

# Find directories with depth 2 in data_dir, this correspond to tissue/condition combinations (e.g. pancreas/tumor)
# Exclude empty directories
tissue_condition_dirs=$(find $data_dir -mindepth 2 -maxdepth 2 -type d -not -empty)

# Init waiting list
waiting_list=""

# Loop over all subtissue/condition directories
for tissue_condition_dir in $tissue_condition_dirs; do
    tissue=$(basename $(dirname $tissue_condition_dir))
    condition=$(basename $tissue_condition_dir)

    echo "$tissue $condition"
    
    if [ $NETWORKS -eq 1 ]; then
        echo 'Networks'

        # create job script
        networks_name="net_${tissue}_${condition}"
        networks_script="$script_dir/$networks_name.sh"

        networks_id=$(fsub \
            -p "$networks_script" \
            -n "$networks_name" \
            -nc "$NETWORKS_NCPUS" \
            -m "$NETWORKS_MEMORY" \
            -e "WGCNA" \
            -q "$NETWORKS_QUEUE" \
            -c "python aggregate_networks.py --inputdir $tissue_condition_dir"
        )

        waiting_list="$waiting_list:$networks_id"
    fi
done

# Directory containing normal and tumor directories with all the tissue networks
tissue_nets_dir="/home/lnemati/pathway_crosstalk/data/networks"
# Find all files in the tissue_nets_dir
networks=$(find $tissue_nets_dir -type f -name *.csv.gz)

# Loop over all tissue networks
for network in $networks; do
    tissue=$(sed 's/.csv.gz//' <<< $(basename $network))
    condition=$(basename $(dirname $network))

    echo "$tissue $condition"
    echo "$network"
    
    if [ $MOTIFS -eq 1 ]; then
        echo 'Motifs'

        # create job script
        motifs_name="motifs_${tissue}_${condition}"
        motifs_script="$script_dir/$motifs_name.sh"

        motifs_id=$(fsub \
            -p "$motifs_script" \
            -n "$motifs_name" \
            -nc "$MOTIFS_NCPUS" \
            -m "$MOTIFS_MEMORY" \
            -e "WGCNA" \
            -q "$MOTIFS_QUEUE" \
            -w "$waiting_list" \
            -c "python motifs.py --inputfile $network"
        )
    fi 
done






#if [ $COEXP -eq 1 ]; then
#    echo 'Coexp'
#
#    # create job script for all tissues
#    coexp_name="coexp"
#    coexp_script="$script_dir/$coexp_name.sh"
#
#    coexp_id=$(fsub \
#        -p "$coexp_script" \
#        -n "$coexp_name" \
#        -nc "$COEXP_NCPUS" \
#        -m "$COEXP_MEMORY" \
#        -e "WGCNA" \
#        -q "$COEXP_QUEUE" \
#        -c "python coexp.py" )
#fi
#
#if [ $MOTIFS -eq 1 ]; then
#    echo 'Motifs'
#
#    # create job script for all tissues
#    motifs_name="motifs"
#    motifs_script="$script_dir/$motifs_name.sh"
#    
#    # if coexp is 1, add coexp to the waiting list 
#    [[ $COEXP -eq 1 ]] && waiting_list="$coexp_id"
#
#    motifs_id=$(fsub \
#        -p "$motifs_script" \
#        -n "$motifs_name" \
#        -nc "$MOTIFS_NCPUS" \
#        -m "$MOTIFS_MEMORY" \
#        -e "WGCNA" \
#        -q "$MOTIFS_QUEUE" \
#        -w "$waiting_list" \
#        -c "python motifs.py" )
#fi
#
#if [ $DENS -eq 1 ]; then
#    echo 'Dens'
#
#    # create job script for all tissues
#    dens_name="dens"
#    dens_script="$script_dir/$dens_name.sh"
#    
#    [[ $COEXP -eq 1 ]] && waiting_list="$coexp_id"
#    [[ $MOTIFS -eq 1 ]] && waiting_list="$motifs_id"
#
#    dens_id=$(fsub \
#        -p "$dens_script" \
#        -n "$dens_name" \
#        -nc "$DENS_NCPUS" \
#        -m "$DENS_MEMORY" \
#        -e "WGCNA" \
#        -q "$DENS_QUEUE" \
#        -w "$waiting_list" \
#        -c "python density.py" )
#fi

echo "Done: pipeline.sh"
