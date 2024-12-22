import pandas as pd
import os
import numpy as np
import gseapy as gp
import argparse
import igraph as ig
from igraph import Graph
from itertools import combinations
from multiprocessing import Pool
import ast
import sys

n_permutations = 5000
ncpus = os.environ['NCPUS']

# Get actual interactions
ccc = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv')
ccc['complex_a'] = ccc['complex_a'].apply(lambda x: ast.literal_eval(x))
ccc['complex_b'] = ccc['complex_b'].apply(lambda x: ast.literal_eval(x))
ccc['all_genes'] = ccc['all_genes'].apply(lambda x: ast.literal_eval(x)).apply(set)
ccc_gene_sets = set(tuple(sorted(gene_set)) for gene_set in ccc['all_genes'])

def find_motifs(
    network,
    adj,
    ccc_gene_sets,
    ):
    
    # Generate pairs of significant interactions
    int_pairs = pd.Series(combinations(pd.Index(network.query('ccc')[['complex1', 'complex2']]), 2))
    int_pairs.index = int_pairs.apply(lambda x: '+'.join(x[0]) + '&' + '+'.join(x[1]))
    int_pairs_genes = int_pairs.apply(lambda x: sorted(set(x[0]).union(x[1])))
    
    # Divide into triplets (2 interactions share a complex) and quadruplets (no shared complexes)
    triplets = int_pairs_genes[int_pairs_genes.apply(len) == 3]
    quadruplets = int_pairs_genes[int_pairs_genes.apply(len) == 4]

    # Count triplets
    g1 = ig.Graph([(0, 1), (1, 2)], directed=False) # 0-1-2 e.g 2 ligands share same receptor, but ligs aren't coexpressed
    g2 = ig.Graph([(0, 1), (1, 2), (2, 0)], directed=False) # clique

    possible_motifs = {
        '3_path': g1,
        '3_clique': g2
    }

    motifs = {}

    for k in possible_motifs.keys():
        motifs[k] = []

    for label, genes in triplets.items():
        subgraph = Graph.Adjacency(adj.loc[genes, genes], mode=ig.ADJ_UNDIRECTED)
        for key, motif in possible_motifs.items():
            if subgraph.isomorphic(motif):
                motifs[key].append(label)
                continue
    
    # Count quadruplets
    g1 = ig.Graph([(0, 3), (1, 2)], directed=False) # no crosstalk
    g2 = ig.Graph([(0, 3), (1, 2), (0, 1)], directed=False) # path or 3-chain
    g3 = ig.Graph([(0, 3), (1, 2), (0, 2), (0, 1)], directed=False) # triangle with exta link
    g4 = ig.Graph([(1, 3), (0, 3), (1, 2), (0, 2)], directed=False) # cycle (4 nodes loop)
    g5 = ig.Graph([(1, 3), (0, 3), (1, 2), (0, 2), (0, 1)], directed=False) # one missing (al pairs are connected except one)
    g6 = ig.Graph([(2, 3), (1, 3), (0, 3), (1, 2), (0, 2), (0, 1)], directed=False) # clique

    possible_motifs = {
        '4_no_crosstalk': g1,
        '4_path': g2,
        '4_triangle_extra': g3,
        '4_cycle' :g4,
        '4_one_missing' :g5,
        '4_clique' :g6,
    }

    for k in possible_motifs.keys():
        motifs[k] = []
        
    for label, genes in quadruplets.items():
        subgraph = Graph.Adjacency(adj.loc[genes, genes], mode=ig.ADJ_UNDIRECTED)
        for key, motif in possible_motifs.items():
            if subgraph.isomorphic(motif):
                motifs[key].append(label)
                continue
    
    # Make series with number of motifs
    counts = pd.Series({k: len(v) for k, v in motifs.items()})
    
    return motifs, counts

def permutation_test(args):
    i, ccc, network = args

    np.random.seed(i)          
    
    network = network.copy()
    ccc = ccc.copy()

    fake_all_genes = ccc.complex_a + ccc.complex_b.sample(frac=1, replace=False).reset_index(drop=True)
    fake_ccc_gene_sets = set(tuple(sorted(gene_set)) for gene_set in fake_all_genes)
    network['ccc'] = False                 
    network['ccc'] = network['all_genes'].apply(lambda genes: tuple(genes) in fake_ccc_gene_sets)

    _, fake_counts = find_motifs(network, adj, fake_ccc_gene_sets)
    return fake_counts  

# Run code 3 times, one for tumor (positive), one for normal (negative), and one for both
for condition in ['tumor', 'normal', 'both']:
    print('Condition:', condition)
    network = pd.read_csv('/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/network_filtered.csv')
    network['all_genes'] = network['all_genes'].apply(lambda x: ast.literal_eval(x))
    
    if condition == 'tumor': # Only keep auroc > 0.5
        network = network.query('auroc > 0.5')
    elif condition == 'normal': # Only keep auroc < 0.5
        network = network.query('auroc < 0.5')

    edges = list(zip(network['complex1'], network['complex2']))
    G = Graph.TupleList(edges, directed=False)

    nodes = G.vs['name']
    adj = pd.DataFrame(G.get_adjacency(), index=nodes, columns=nodes)

    # Find motifs
    motifs, counts = find_motifs(network, adj, ccc_gene_sets)

    # Save results
    outdir = '/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/motifs/' + condition
    os.makedirs(outdir, exist_ok=True)
    # set counts columns to motif and number
    counts = counts.reset_index()
    counts.columns = ['motif', 'number']
    counts.to_csv(os.path.join(outdir, 'counts.csv'), index=False)
    # Motifs is a dictionary of lists of interactions that form each motif,
    # make a dataframe with one row per interaction and one column for motif type
    data = [(motif_type, interaction) for motif_type, interactions in motifs.items() for interaction in interactions]
    df = pd.DataFrame(data, columns=['Type', 'Interaction'])
    df.to_csv(os.path.join(outdir, 'motifs.csv'), index=False)
    print('Actual motifs detected:', counts)

    # Permutation testing using multiprocessing
    print('Running permutations')
    with Pool(int(ncpus)) as p:
        fake_counts = pd.DataFrame(p.map(permutation_test, [(i, ccc, network) for i in range(n_permutations)]))

    print('Saving results')
    fake_counts.to_csv(os.path.join(outdir, 'motifs_permutations.csv'), index=False)

print('Done: motifs.py')
