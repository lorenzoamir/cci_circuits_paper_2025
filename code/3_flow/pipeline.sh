#!/bin/bash
#PBS -N wgcna_analysis
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

FLOW=1
GENERATE=0 # Generate bootstrap interactions
BOOTSTRAP=0 # Get same_module for bootstrap interactions 
AGGREGATE=0 # Aggregate bootstrap results
SANKEY=0

FLOW_QUEUE='q02anacreon'
GENERATE_QUEUE='q02anacreon'
BOOTSTRAP_QUEUE='q02anacreon'
AGGREGATE_QUEUE='q02anacreon'
SANKEY_QUEUE='q02anacreon'

FLOW_NCPUS=8
GENERATE_NCPUS=8
BOOTSTRAP_NCPUS=2
AGGREGATE_NCPUS=8
SANKEY_NCPUS=8

FLOW_MEMORY=8gb
GENERATE_MEMORY=20gb
BOOTSTRAP_MEMORY=8gb
AGGREGATE_MEMORY=16gb
SANKEY_MEMORY=8gb

cd /home/lnemati/pathway_crosstalk/code/3_flow
script_dir="/home/lnemati/pathway_crosstalk/code/3_flow/scripts"
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
if [ $FLOW -eq 1 ]; then
    echo 'Flow'

    # create job script for all tissues
    flow_name="flow"
    flow_script="$script_dir/$flow_name.sh"

    echo "$flow_script"
    
    flow_id=$(fsub \
        -p "$flow_script" \
        -n "$flow_name" \
        -nc "$FLOW_NCPUS" \
        -m "$FLOW_MEMORY" \
        -e "WGCNA" \
        -q "$FLOW_QUEUE" \
        -c "python flow.py --tissues \"$tissues_string\"")

fi

if [ $GENERATE -eq 1 ]; then
    echo 'Generate'

    # create job script for all tissues
    generate_name="generate"
    generate_script="$script_dir/$generate_name.sh"

    generate_id=$(fsub \
        -p "$generate_script" \
        -n "$generate_name" \
        -nc "$GENERATE_NCPUS" \
        -m "$GENERATE_MEMORY" \
        -e "WGCNA" \
        -q "$GENERATE_QUEUE" \
        -c "python generate_bootstrap.py --tissues \"$tissues_string\"")
fi

# Generate list of wgcna files searching for wgcna_*.p
mapfile -t wgcna_files < <(find "$data_dir" -name "wgcna_*.p")

for wgcnafile in "${wgcna_files[@]}"; do
    if [ $BOOTSTRAP -eq 1 ]; then
        # Cut 'wgcna_' and '.p' from the file name
        tissue_name=$(basename "$wgcnafile" | sed 's/wgcna_//g' | sed 's/\.p//g')

        # create job script for each tissue
        bootstrap_name="btp_$tissue_name"
        bootstrap_script="$script_dir/$bootstrap_name.sh"
        
        # Wait for generate to finish
        waiting_list=""
        [ $GENERATE -eq 1 ] && waiting_list="$waiting_list:$generate_id"

        bootstrap_id=$(fsub \
            -p "$bootstrap_script" \
            -n "$bootstrap_name" \
            -nc "$BOOTSTRAP_NCPUS" \
            -m "$BOOTSTRAP_MEMORY" \
            -e "WGCNA" \
            -q "$BOOTSTRAP_QUEUE" \
            -w "$waiting_list" \
            -c "python bootstrap_interactions.py --input $wgcnafile")
    fi
done

# AGGREAGATE takes all tissues like FLOW
if [ $AGGREGATE -eq 1 ]; then
    echo 'Aggregate'

    # create job script for all tissues
    aggregate_name="aggregate"
    aggregate_script="$script_dir/$aggregate_name.sh"

    # Wait for generate and bootstrap
    waiting_list=""
    [ $GENERATE -eq 1 ] && waiting_list="$waiting_list:$generate_id"
    [ $BOOTSTRAP -eq 1 ] && waiting_list="$waiting_list:$bootstrap_id"

    aggregate_id=$(fsub \
        -p "$aggregate_script" \
        -n "$aggregate_name" \
        -nc "$AGGREGATE_NCPUS" \
        -m "$AGGREGATE_MEMORY" \
        -e "WGCNA" \
        -q "$AGGREGATE_QUEUE" \
        -c "python aggregate_bootstrap.py --tissues \"$tissues_string\"")
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

echo "Done: pipeline.sh"
