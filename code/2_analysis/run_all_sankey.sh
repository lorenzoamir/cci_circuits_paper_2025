#!/bin/bash
#PBS -N sankey
#PBS -l select=1:ncpus=2:ngpus=0:mem=8gb
#PBS -q q02anacreon

source /projects/bioinformatics/snsutils/snsutils.sh

# Define tumor paths
tumor_paths=(
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/uterus/tumor/uterine_carcinosarcoma"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/thyroid/tumor/thyroid_carcinoma"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/kidney/tumor/kidney_papillary_cell_carcinoma"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/kidney/tumor/kidney_chromophobe"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/kidney/tumor/kidney_clear_cell_carcinoma"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/adrenal_gland/tumor/adrenocortical_cancer"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/liver/tumor/liver_hepatocellular_carcinoma"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/prostate/tumor/prostate_adenocarcinoma"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/lung/tumor/lung_squamous_cell_carcinoma"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/lung/tumor/lung_adenocarcinoma"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/breast/tumor/breast_invasive_carcinoma"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/ovary/tumor/ovarian_serous_cystadenocarcinoma"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/pancreas/tumor/pancreatic_adenocarcinoma"
    "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/stomach/tumor/stomach_adenocarcinoma"
)

# Activate conda environment
cd /home/lnemati/pathway_crosstalk/code/2_analysis
eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA

# Iterate through tumor paths
for tumor_path in "${tumor_paths[@]}"; do
    # Derive tumor file path
    tumor_file="$tumor_path/interactions/ccc.csv"
    
    # Derive normal file path
    # Replace "/tumor/" with "/normal/", and replace the last segment with the tissue name
    tissue_name=$(basename "$(dirname "$(dirname "$tumor_path")")")
    # Replace /tumor/ with /normal/
    normal_dir=$(echo "$tumor_file" | sed "s|/tumor/.*|/normal/|")
    # Normal tissue name is the only directory contained in normal_dir
    normal_tissue_name=$(basename $(find "$normal_dir" -maxdepth 1 -mindepth 1 -type d))
    # Add tissue_name/interactions/ccc.csv
    normal_file="$normal_dir$normal_tissue_name/interactions/ccc.csv"
    
    # Extract tissue name for logging
    
    # Print the details
    echo "$tissue_name"
    echo "python sankey.py --tumor $tumor_file --normal $normal_file"
    
    id=$(fsub \
        -p "/home/lnemati/pathway_crosstalk/code/2_analysis/scripts/sankey_$tissue_name.sh" \
        -n "sankey_$tissue_name" \
        -q q02gaia \
        -nc 1 \
        -m '8gb' \
        -e 'WGCNA' \
        -c "python sankey.py --tumor $tumor_file --normal $normal_file")

    echo ""
done

exit 0

