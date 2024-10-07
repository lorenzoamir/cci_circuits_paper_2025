import numpy as np
import pandas as pd
import sys
import os
import argparse
import gseapy as gp
from itertools import combinations

# read filename from command line
parser = argparse.ArgumentParser(description='Enrichment analysis of modules')
parser.add_argument('-i', '--input', type=str, help='path tocsv file containing modules', required=True)
args = parser.parse_args()

# Output path
output_path = os.path.dirname(args.input)

# Read modules
modules = pd.read_csv(args.input, index_col=0)
# Make into a series
modules = modules["module"]
all_genes = modules.index.tolist()

reactome = '/home/lnemati/resources/reactome/ReactomePathways.gmt'


# ----- Enrichment analysis ------

print("Enrichment analysis")
enrichment = pd.DataFrame()

for module in modules.unique():
    module_genes = modules[modules == module].index.tolist()
    # if less than 5 genes in module, skip
    if len(module_genes) < 5:
        continue
    try:
        enr = gp.enrichr(
            gene_list=module_genes,
            #gene_sets="Reactome_2022",
            gene_sets=reactome,
            background=all_genes,
            outdir=None,
            )
    except Exception as e:
        print("Error in enrichment analysis for module: ", module, file=sys.stderr)
        print("Error in enrichment analysis for module: ", module, file=sys.stdout)
        print()
        print("Exception:", file=sys.stdout)
        print("Exception:", file=sys.stderr)
        print(e, file=sys.stdout)
        print(e, file=sys.stderr)
        print("Module size: ", len(module_genes), file=sys.stderr)
        print("Module size: ", len(module_genes), file=sys.stdout)
        continue
    enr = enr.results
    # Only keep significant enrichments
    enr = enr[enr["Adjusted P-value"] < 0.05]
    # Add module column and concatenate
    enr["module"] = module
    enrichment = pd.concat([enrichment, enr])

# Save to csv
enrichment.to_csv(os.path.join(output_path, 'enrichment.csv'), index=False)

print("Done: enrichment.py")
