import numpy as np
import sys
import PyWGCNA
import pandas as pd
import os
import argparse
import gseapy as gp
from itertools import combinations

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

#reactome = '/home/lnemati/pathway_crosstalk/results/reactome/reactome_connected_pathways.gmt'
reactome = '/home/lnemati/resources/reactome/ReactomePathways.gmt'

# ----- Enrichment analysis ------

print("Enrichment analysis")
enrichment = pd.DataFrame()

for module in WGCNA.datExpr.var["moduleLabels"].unique():
    module_genes = WGCNA.datExpr.var[WGCNA.datExpr.var["moduleLabels"] == module].index.tolist()
    # if less than 5 genes in module, skip
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
        print("Module size: ", len(module_genes), file=sys.stderr)
        print("Module size: ", len(module_genes), file=sys.stdout)
        print("Error: ", e, file=sys.stderr)
        print("Error: ", e, file=sys.stdout)
        continue
    enr = enr.results
    # Only keep significant enrichments
    # enr = enr[enr["Adjusted P-value"] < 0.05]
    # Add module column and concatenate
    enr["module"] = module
    enrichment = pd.concat([enrichment, enr])

print("Saving enrichment results to: ", os.path.join(output_path, "module_enrichment.csv.gz"))
enrichment.to_csv(os.path.join(output_path, "module_enrichment.csv.gz"), index=False)

#all_terms = enrichment["Term"].unique()
## ----- Occurrences -----
#
## Count occurrences of terms in modules
#terms_counts = pd.DataFrame(0, index=all_terms, columns=WGCNA.datExpr.var["moduleLabels"].unique())
#
#for module, df in enrichment.groupby("module"):
#    for term in df["Term"]:
#        terms_counts.loc[term, module] += 1
#
#print("Saving terms occurrences to: ", os.path.join(output_path, "pathways_occurrences.csv.gz"))
#terms_counts.to_csv(os.path.join(output_path, "pathways_occurrences.csv.gz"))
#
## ---- Co-occurrences -----
#
## Count co-occurrences of terms in the same module
#cooccurrences = pd.DataFrame(0, index=all_terms, columns=all_terms)
#
#for module in WGCNA.datExpr.var["moduleLabels"].unique():
#    oc = terms_counts[terms_counts[module] > 0]
#    pairs = list(combinations(oc.index, 2))
#    for pair in pairs:
#        cooccurrences.loc[pair[0], pair[1]] += 1
#        cooccurrences.loc[pair[1], pair[0]] += 1
#
#print("Saving terms co-occurrences to: ", os.path.join(output_path, "pathways_cooccurrences.csv.gz"))
#cooccurrences.to_csv(os.path.join(output_path, "pathways_cooccurrences.csv.gz"))
#
## Count co-occurrences of terms in the same module
#terms_counts = pd.DataFrame(0, index=all_terms, columns=all_terms)
#
#for module, df in enrichment.groupby("module"):
#    # generate pairwise combinations of terms, no repetitions
#    pairs = list(combinations(df["Term"], 2))
#    for pair in pairs:
#        terms_counts.loc[pair[0], pair[1]] += 1
#        terms_counts.loc[pair[1], pair[0]] += 1
#
#print("Saving terms co-occurrences to: ", os.path.join(output_path, "pathways_cooccurrences.csv.gz"))
#terms_counts.to_csv(os.path.join(output_path, "pathways_cooccurrences.csv.gz"))

print("Done: enrichment.py")
