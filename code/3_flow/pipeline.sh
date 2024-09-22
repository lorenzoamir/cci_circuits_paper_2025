#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

FLOW=1
SANKEY=0

FLOW_QUEUE='q02anacreon'
SANKEY_QUEUE='q02anacreon'

FLOW_NCPUS=8
SANKEY_NCPUS=8

FLOW_MEMORY=8gb
SANKEY_MEMORY=8gb

cd /home/lnemati/pathway_crosstalk/code/3_flow
script_dir="/home/lnemati/pathway_crosstalk/code/3_flow/scripts"
if [ ! -d "$script_dir" ]; then
    mkdir "$script_dir"
fi

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

# find all corresponding_normal_wgcna.txt files in the data directory
data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/

# Find all tissue directories
mapfile -t tissuedirs < <(find "$data_dir" -maxdepth 1 -mindepth 1 -type d)

tumors=()
normals=()

for tissuedir in "${tissuedirs[@]}"; do
    # Inside each tissue directory, there are tumor and normal directories
    tumordir="$tissuedir/tumor"
    normaldir="$tissuedir/normal"

    # Check if both tumor and normal directories contain one and only one directory each
    tumorcount=$(find "$tumordir" -maxdepth 1 -mindepth 1 -type d | wc -l)
    normalcount=$(find "$normaldir" -maxdepth 1 -mindepth 1 -type d | wc -l)

    # Also check that none of the directories they contain are empty
    t_empty=$(find "$tumordir" -maxdepth 1 -mindepth 1 -type d -empty | wc -l)
    n_empty=$(find "$normaldir" -maxdepth 1 -mindepth 1 -type d -empty | wc -l)

    if [ "$tumorcount" -eq 1 ] && [ "$normalcount" -eq 1 ] && [ "$t_empty" -eq 0 ] && [ "$n_empty" -eq 0 ]; then
        tumordir=$(find "$tumordir" -maxdepth 1 -mindepth 1 -type d)
        normaldir=$(find "$normaldir" -maxdepth 1 -mindepth 1 -type d)

        # Save the tumor and normal paths
        tumors+=("$tumordir")
        normals+=("$normaldir")
        echo "$tumordir"
        echo "$normaldir"
        echo ""
    fi
done

# FLOW takes all tissues at once
if [ $FLOW -eq 1 ]; then
    echo 'Flow'

    # create job script for all tissues
    flow_name="flow"
    flow_script="$script_dir/$flow_name.sh"

    echo "$flow_script"
    
    tumors_str=\"${tumors[@]}\"
    normals_str=\"${normals[@]}\"
    
    echo tumors_str
    echo normals_str

    flow_id=$(fsub \
        -p "$flow_script" \
        -n "$flow_name" \
        -nc "$FLOW_NCPUS" \
        -m "$FLOW_MEMORY" \
        -e "WGCNA" \
        -q "$FLOW_QUEUE" \
        -c "python flow.py --tumors $tumors_str --normals $normals_str")

fi

# Run pipeline for each tissue
for i in "${!tumors[@]}"; do
    tumordir="${tumors[$i]}"
    normaldir="${normals[$i]}"

    # Run pipeline for each tissue
    if [ $SANKEY -eq 1 ]; then
        echo 'Sankey'

        # create job script for each tumor
        sankey_name="sankey_$name"
        sankey_script="$tumordir/scripts/sankey.sh"

        sankey_id=$(fsub \
            -p "$sankey_script" \
            -n "$sankey_name" \
            -nc "$SANKEY_NCPUS" \
            -m "$SANKEY_MEMORY" \
            -e "WGCNA" \
            -q "$SANKEY_QUEUE" \
            -c "python sankey.py --tumordir $tumordir --normaldir $normaldir")
    fi
done
