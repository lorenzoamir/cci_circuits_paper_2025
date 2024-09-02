#!/bin/bash
#PBS -N coevolution_pipeline
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source "/projects/bioinformatics/snsutils/snsutils.sh"

MUTATIONS=0
SEPARATE=0
JACCARD=1
COOCCURRENCE=0
DIRECTIONALITY=0

MUTATIONS_QUEUE="q02anacreon"
SEPARATE_QUEUE="q02anacreon"
JACCARD_QUEUE="q02anacreon"
COOCCURRENCE_QUEUE="q02anacreon"
DIRECTIONALITY_QUEUE="q02anacreon"

MUTATIONS_NCPUS=16
SEPARATE_NCPUS=8
JACCARD_NCPUS=32
COOCCURRENCE_NCPUS=8
DIRECTIONALITY_NCPUS=8

MUTATIONS_MEMORY=16gb
SEPARATE_MEMORY=8gb
JACCARD_MEMORY=32gb
COOCCURRENCE_MEMORY=16gb
DIRECTIONALITY_MEMORY=16gb

lr_resource="/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv"

cd "/home/lnemati/pathway_crosstalk/code/1_coevolution"

# make scripts_dir
scripts_dir="/home/lnemati/pathway_crosstalk/code/1_coevolution/scripts"
mkdir -p "$scripts_dir"

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

# ----- MUTATIONS -----

if [ $MUTATIONS -eq 1 ]; then
    echo "MUTATIONS"

    mutations_name=mutations
    mutations_script="$scripts_dir"/"$mutations_name".sh

    mutations_id=$(fsub \
        -p "$mutations_script" \
        -n "$mutations_name" \
        -nc "$MUTATIONS_NCPUS" \
        -m "$MUTATIONS_MEMORY" \
        -q "$MUTATIONS_QUEUE" \
        -e "WGCNA" \
        -c "python get_mutations.py")
fi

if [ $SEPARATE -eq 1 ]; then
    echo "SEPARATE"

    separate_name=separate
    separate_script="$scripts_dir"/"$separate_name".sh

    # Wait for mutations job to finish
    waiting_list=""
    [ $MUTATIONS -eq 1 ] && waiting_list="$waiting_list:$mutations_id"

    separate_id=$(fsub \
        -p "$separate_script" \
        -n "$separate_name" \
        -nc "$SEPARATE_NCPUS" \
        -m "$SEPARATE_MEMORY" \
        -q "$SEPARATE_QUEUE" \
        -e "WGCNA" \
        -w "$waiting_list" \
        -c "python separate.py")
fi

if [ $JACCARD -eq 1 ]; then
    echo "JACCARD"

    jaccard_name=jacc
    jaccard_script="$scripts_dir"/"$jaccard_name".sh

    inputfile="/home/lnemati/pathway_crosstalk/data/tumor_coev/all_cancers/Patient/mutations_df.csv"

    # Wait for separate job to finish
    waiting_list=""
    [ $MUTATIONS -eq 1 ] && waiting_list="$waiting_list:$mutations_id"
    [ $SEPARATE -eq 1 ] && waiting_list="$waiting_list:$separate_id"

    jaccard_id=$(fsub \
        -p "$jaccard_script" \
        -n "$jaccard_name" \
        -nc "$JACCARD_NCPUS" \
        -m "$JACCARD_MEMORY" \
        -q "$JACCARD_QUEUE" \
        -e "WGCNA" \
        -w "$waiting_list" \
        -c "python jaccard.py --inputfile $inputfile")
fi

datadir="/home/lnemati/pathway_crosstalk/data/tumor_coev"

# Find all cancertype (tissue) directories in data directory
# Structure is datadir/cancertype/mutationtype/mutations_df.csv
# e.g. /home/lnemati/pathway_crosstalk/data/tumor_coev/all_cancers/Private/mutations_df.csv

tissuedirs=($(find "$datadir" -mindepth 1 -maxdepth 1 -type d))

for tissuedir in "${tissuedirs[@]}"; do
    tissue=$(basename "$tissuedir")
    echo "$tissue"

    patient="$tissuedir/Patient/mutations_df.csv"
    trunk="$tissuedir/Trunk/mutations_df.csv"
    branch="$tissuedir/Branch/mutations_df.csv"
    private="$tissuedir/Private/mutations_df.csv"

    if [ $COOCCURRENCE -eq 1 ]; then
        echo "COOCCURRENCE:"

        cooccurrence_name=cooc_"$tissue"
        cooccurrence_script="$scripts_dir"/"$cooccurrence_name".sh

        # Wait for separate job to finish
        waiting_list=""
        [ $MUTATIONS -eq 1 ] && waiting_list="$waiting_list:$mutations_id"
        [ $SEPARATE -eq 1 ] && waiting_list="$waiting_list:$separate_id"

        cooccurrence_id=$(fsub \
            -p "$cooccurrence_script" \
            -n "$cooccurrence_name" \
            -nc "$COOCCURRENCE_NCPUS" \
            -m "$COOCCURRENCE_MEMORY" \
            -q "$COOCCURRENCE_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python cooccurrence.py --inputfile $patient")
    fi

    if [ $DIRECTIONALITY -eq 1 ]; then
        echo "DIRECTIONALITY:"

        directional_name=dir_"$tissue"
        directional_script="$scripts_dir"/"$directional_name".sh

        # Wait for separate job to finish
        waiting_list=""
        [ $MUTATIONS -eq 1 ] && waiting_list="$waiting_list:$mutations_id"
        [ $SEPARATE -eq 1 ] && waiting_list="$waiting_list:$separate_id"
        [ $COOCCURRENCE -eq 1 ] && waiting_list="$waiting_list:$cooccurrence_id"

        directional_id=$(fsub \
            -p "$directional_script" \
            -n "$directional_name" \
            -nc "$DIRECTIONALITY_NCPUS" \
            -m "$DIRECTIONALITY_MEMORY" \
            -q "$DIRECTIONALITY_QUEUE" \
            -e "WGCNA" \
            -w "$waiting_list" \
            -c "python directionality.py --mutationsdir $tissuedir --significant $tissuedir/Patient/cooccurrences_filtered.csv")
    fi

echo ""
done

echo "Done: pipeline.sh"

## ----- JACCARD -----
#
#if [ $JACCARD -eq 1 ]; then
#    echo "JACCARD"
#
#    jaccard_name=jaccard
#    jaccard_script="$scripts_dir"/"$jaccard_name".sh
#
#    jaccard_id=$(fsub \
#        -p "$jaccard_script" \
#        -n "$jaccard_name" \
#        -nc "$JACCARD_NCPUS" \
#        -m "$JACCARD_MEMORY" \
#        -q "$JACCARD_QUEUE" \
#        -e "WGCNA" \
#        -c "python coevolution_jaccard.py")
#
#fi
#
## ----- COOCCURRENCE -----
#if [ $COOCCURRENCE -eq 1 ]; then
#    echo "COOCCURRENCE"
#
#    cooccurrence_name=cooccurrence
#    cooccurrence_script="$scripts_dir"/"$cooccurrence_name".sh
#
#    cooccurrence_id=$(fsub \
#        -p "$cooccurrence_script" \
#        -n "$cooccurrence_name" \
#        -nc "$COOCCURRENCE_NCPUS" \
#        -m "$COOCCURRENCE_MEMORY" \
#        -q "$COOCCURRENCE_QUEUE" \
#        -e "WGCNA" \
#        -c "python cooccurrence.py")
#fi
#
## ----- DIRECTIONALITY -----
#if [ $DIRECTIONALITY -eq 1 ]; then
#    echo "DIRECTIONALITY"
#
#    directional_name=directionality
#    directional_script="$scripts_dir"/"$directional_name".sh
#
#    directional_id=$(fsub \
#        -p "$directional_script" \
#        -n "$directional_name" \
#        -nc "$DIRECTIONALITY_NCPUS" \
#        -m "$DIRECTIONALITY_MEMORY" \
#        -q "$DIRECTIONALITY_QUEUE" \
#        -e "WGCNA" \
#        -c "python directionality.py")
#fi
#
#exit 0
#
#
