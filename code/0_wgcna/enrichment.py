import numpy as np
import sys
import PyWGCNA
import pandas as pd
import os
import argparse
import gseapy as gp

# read filename from command line
parser = argparse.ArgumentParser(description='Enrichment analysis of WGCNA modules.')
parser.add_argument('-i', '--input', type=str, help='path to input WGCNA object.')
args = parser.parse_args()

# Output path
output_path = os.path.dirname(args.input)

# Load WGCNA object
print("Loading WGCNA object: ", args.input)
WGCNA = PyWGCNA.readWGCNA(args.input)
all_genes = WGCNA.datExpr.var.index.tolist()
# ----- Enrichment analysis ------

print("Enrichment analysis")
enrichment = pd.DataFrame()

for module in WGCNA.datExpr.var["moduleLabels"].unique():
    module_genes = WGCNA.datExpr.var[WGCNA.datExpr.var["moduleLabels"] == module].index.tolist()
    enr = gp.enrichr(
        gene_list=module_genes,
        gene_sets="Reactome_2022",
        background=all_genes,
        outdir=None,
        )
    enr = enr.results
    # Only keep significant enrichments
    enr = enr[enr["Adjusted P-value"] < 0.05]
    # Add module column and concatenate
    enr["module"] = module
    enrichment = pd.concat([enrichment, enr])

print("Saving enrichment results to: ", os.path.join(output_path, "module_enrichment.csv.gz"))
enrichment.to_csv(os.path.join(output_path, "module_enrichment.csv.gz"), index=False)

print("Done: enrichment.py")
