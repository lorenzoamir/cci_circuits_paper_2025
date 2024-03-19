#!/bin/bash
#PBS -N wgcna_master
#PBS -l select=1:ncpus=32:mem=64gb
#PBS -q q02anacreon

SEPARATE=0

DESEQ=0
DESEQ_QUEUE="q02gaia"
DESEQ_NCPUS=8
DESEQ_MEMORY=16gb

WGCNA=0
WGCNA_QUEUE="q02anacreon"
WGCNA_NCPUS=8
WGCNA_MEMORY=20gb

CCC=0
CCC_QUEUE="q02anacreon"
CCC_NCPUS=1
CCC_MEMORY=20gb

ENRICHMENT=1
ENRICHMENT_QUEUE="q02gaia"
ENRICHMENT_NCPUS=2
ENRICHMENT_MEMORY=8gb

# get ncpus from NCPUS environment variable
ncpus=8
memory=20gb
queue="q02anacreon"

lr_resource="/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/consensus.csv"

cd "/home/lnemati/pathway_crosstalk/code/0_WGCNA_CCC"

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"

# ----- SEPARATE -----
if [ $SEPARATE -eq 1 ]; then
    conda activate WGCNA
    # separate data according to tissue 
    echo "Separating data according to tissue"
    echo ""
    python /home/lnemati/pathway_crosstalk/code/0_WGCNA_CCC/SeparateTissuesConditions.py --inputfile "/projects/bioinformatics/DB/Xena/TCGA_GTEX/TCGA_GTEX.h5ad" --outputdir "/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/"

    echo "Done"
    echo ""
fi

# ----- Analysis -----

# find all .h5ad files in the data directory and its subdirectories
data_dir=/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/
mapfile -t files < <(find "$data_dir" -name "*.h5ad")

eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"
conda activate WGCNA
 
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

        echo "Creating job script for $deseq_script"
        touch "$deseq_script"

        echo "#!/bin/bash" > "$deseq_script"
        echo "#PBS -l select=1:ncpus=$DESEQ_NCPUS:mem=$DESEQ_MEMORY" >> "$deseq_script"
        echo "#PBS -q $DESEQ_QUEUE" >> "$deseq_script"
        echo "#PBS -N $deseq_name" >> "$deseq_script"
        echo "" >> "$deseq_script"
        echo 'eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"' >> "$deseq_script"
        echo "conda activate WGCNA" >> "$deseq_script"
        echo "python /home/lnemati/pathway_crosstalk/code/0_WGCNA_CCC/deseq_norm.py --input "${file} --ncpus $DESEQ_NCPUS"" >> "$deseq_script"
        echo "exit 0" >> "$deseq_script"

        # submit job
        echo "Submitting job for $deseq_script"
        deseq_id=$(qsub "$deseq_script")

        # add deseq_id to waiting list
        waiting_list="$waiting_list:$deseq_id"
        echo "Waiting list: $waiting_list"
        echo ""
    fi

    # ----- WGCNA -----
    if [ $WGCNA -eq 1 ]; then
        echo "WGCNA"

        # Create job script
        wgcna_name=wgcna_"$job_name"
        wgcna_script=$(dirname "$file")/scripts/$wgcna_name.sh
        echo "Creating job script for $wgcna_script"
        touch "$wgcna_script"

        echo "#!/bin/bash" > "$wgcna_script"
        echo "#PBS -l select=1:ncpus=$WGCNA_NCPUS:mem=$WGCNA_MEMORY" >> "$wgcna_script"
        echo "#PBS -q $WGCNA_QUEUE" >> "$wgcna_script"
        echo "#PBS -N $wgcna_name" >> "$wgcna_script"
        echo "" >> "$wgcna_script"
        echo 'eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"' >> "$wgcna_script"
        echo "conda activate WGCNA" >> "$wgcna_script"
        echo "python /home/lnemati/pathway_crosstalk/code/0_WGCNA_CCC/wgcna.py --input "${file}" --ncpus ${ncpus} --lr_resource "${lr_resource}"" >> "${wgcna_script}" 
        echo "exit 0" >> "$wgcna_script"

        # submit job
        echo "Submitting job for $wgcna_script"
        if [ -n "$waiting_list" ]; then
            wgcna_id=$(qsub -W depend=afterok$waiting_list "$wgcna_script")
        else
            wgcna_id=$(qsub "$wgcna_script")
        fi
        # add wgcna_id to waiting list
        waiting_list="$waiting_list:$wgcna_id"
        echo "Waiting list: $waiting_list"
        echo ""
    fi

    # ----- CCC -----
    if [ $CCC -eq 1 ]; then
        echo "CCC"

        # Create job script
        ccc_name=ccc_"$job_name"
        ccc_script=$(dirname "$file")/scripts/$ccc_name.sh
        echo "Creating job script for $ccc_script"
        touch "$ccc_script"

        echo "#!/bin/bash" > "$ccc_script"
        echo "#PBS -l select=1:ncpus=$CCC_NCPUS:mem=$CCC_MEMORY" >> "$ccc_script"
        echo "#PBS -q $CCC_QUEUE" >> "$ccc_script"
        echo "#PBS -N $ccc_name" >> "$ccc_script"
        echo "" >> "$ccc_script"
        echo 'eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"' >> "$ccc_script"
        echo "conda activate c2c" >> "$ccc_script"
        echo "python /home/lnemati/pathway_crosstalk/code/0_WGCNA_CCC/ccc.py --input "${file}"" >> "${ccc_script}"
        echo "exit 0" >> "$ccc_script"

        # submit job
        echo "Submitting job for $ccc_script"
        if [ -n "$waiting_list" ]; then
            qsub -W depend=afterok$waiting_list "$ccc_script"
        else
            qsub "$ccc_script"
        fi
        echo ""
    fi

    # ----- Enrichment analysis -----
    if [ $ENRICHMENT -eq 1 ]; then
        echo "ENRICHMENT"
        # input file is the wgcna output file
        wgcnafile=$(dirname "$file")/WGCNA_"$job_name".p
        echo "Input file for enrichment: $wgcnafile"

        # Create job script
        enrichment_name=enrichment_"$job_name"
        enrichment_script=$(dirname "$file")/scripts/$enrichment_name.sh
        echo "Creating job script for $enrichment_script"
        touch "$enrichment_script"

        echo "#!/bin/bash" > "$enrichment_script"
        echo "#PBS -l select=1:ncpus=$ENRICHMENT_NCPUS:mem=$ENRICHMENT_MEMORY" >> "$enrichment_script"
        echo "#PBS -q $ENRICHMENT_QUEUE" >> "$enrichment_script"
        echo "#PBS -N $enrichment_name" >> "$enrichment_script"
        echo "" >> "$enrichment_script"
        echo 'eval "$(/cluster/shared/software/miniconda3/bin/conda shell.bash hook)"' >> "$enrichment_script"
        echo "conda activate WGCNA" >> "$enrichment_script" 
        echo "python /home/lnemati/pathway_crosstalk/code/0_WGCNA_CCC/enrichment.py --input "${wgcnafile}"" >> "${enrichment_script}" 
        echo "exit 0" >> "$enrichment_script"

        # submit job
        echo "Submitting job for $enrichment_script"
        if [ -n "$waiting_list" ]; then
            echo "qsub -W depend=afterok$waiting_list $enrichment_script"
            qsub -W depend=afterok$waiting_list "$enrichment_script"
        else
            qsub "$enrichment_script"
        fi
        echo ""
    fi
done

echo "Done: all jobs submitted"

exit 0

