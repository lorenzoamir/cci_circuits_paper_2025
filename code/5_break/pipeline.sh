#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

BREAK=1

BREAK_QUEUE='q02gaia'

BREAK_NCPUS=32

BREAK_MEMORY=100gb

cd /home/lnemati/pathway_crosstalk/code/5_break
script_dir="/home/lnemati/pathway_crosstalk/code/5_break/scripts"
if [ ! -d "$script_dir" ]; then
    mkdir "$script_dir"
fi

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/

# Find all tissue directories
mapfile -t tissuedirs < <(find "$data_dir" -maxdepth 1 -mindepth 1 -type d)

tumors=() # Tumor types with only one normal tissue associated
normals=() # Normal tissues associated with tumor types
tissues_with_both=() # Tissues at least one tumor and one normal tissue

for tissuedir in "${tissuedirs[@]}"; do
    # Inside each tissue directory, there are tumor and normal directories
    tumordir="$tissuedir/tumor"
    normaldir="$tissuedir/normal"

    # Check if normal directories contain one and only one directory each
    tumorcount=$(find "$tumordir" -maxdepth 1 -mindepth 1 -type d | wc -l)
    normalcount=$(find "$normaldir" -maxdepth 1 -mindepth 1 -type d | wc -l)

    # Also check that none of the directories they contain are empty
    t_empty=$(find "$tumordir" -maxdepth 1 -mindepth 1 -type d -empty | wc -l)
    n_empty=$(find "$normaldir" -maxdepth 1 -mindepth 1 -type d -empty | wc -l)
    
    # Check if at least one tumor and one normal tissue are present
    if [ "$tumorcount" -gt 0 ] && [ "$normalcount" -gt 0 ] && [ "$t_empty" -eq 0 ] && [ "$n_empty" -eq 0 ]; then
        tissues_with_both+=("$tissuedir")
        echo "$tissuedir"
        
        # Check if there is only one tumor and one normal tissue
        if [ "$normalcount" -eq 1 ] && [ "$tumorcount" -eq 1 ]; then
            tumordirs=$(find "$tumordir" -maxdepth 1 -mindepth 1 -type d)
            normaldir=$(find "$normaldir" -maxdepth 1 -mindepth 1 -type d)

            # Save the tumor and normal paths pairs (multiple tumors can have the same normal)
            for tumor in $tumordirs; do
                tumors+=("$tumor")
                normals+=("$normaldir")
                echo "$tumor"
                echo "$normaldir"
            done
        fi
    echo ""
    fi
done

# Convert list to string for passing as argument
# The entries are separated by a space, add double quotes to keep them together
tissues_string=$(printf "\"%s\" " "${tissues_with_both[@]}")
echo "$tissues_string"

# FLOW takes all tissues at once
if [ $BREAK -eq 1 ]; then
    echo 'Break'

    # create job script for all tissues
    break_name="bkp"
    break_script="$script_dir/$break_name.sh"

    break_id=$(fsub \
        -p "$break_script" \
        -n "$break_name" \
        -nc "$BREAK_NCPUS" \
        -m "$BREAK_MEMORY" \
        -e "WGCNA" \
        -q "$BREAK_QUEUE" \
        -c "python breakpoints.py --tissues \"$tissues_string\"")
fi

echo "Done: pipeline.sh"
