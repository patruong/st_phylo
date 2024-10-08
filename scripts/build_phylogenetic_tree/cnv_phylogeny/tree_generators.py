import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import pdist

import random
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

def bifurcating_tree(taxa):
    """Create a bifurcating tree with the given taxa."""
    if len(taxa) < 2:
        raise ValueError("At least two taxa are required to create a bifurcating tree.")

    def create_clade(subtaxa):
        if len(subtaxa) == 1:
            return Clade(name=subtaxa[0], branch_length=random.random())
        elif len(subtaxa) == 2:
            clade = Clade(branch_length=random.random())
            clade.clades = [
                Clade(name=subtaxa[0], branch_length=random.random()),
                Clade(name=subtaxa[1], branch_length=random.random())
            ]
            return clade
        else:
            clade = Clade(branch_length=random.random())
            mid = len(subtaxa) // 2
            clade.clades = [create_clade(subtaxa[:mid]), create_clade(subtaxa[mid:])]
            return clade

    tree = Tree(create_clade(taxa))
    return tree

def random_tree(taxa):
    """Create a random binary tree with the given taxa."""
    tree = Tree(Clade())
    clades = [Clade(name=taxon, branch_length=random.random()) for taxon in taxa]
    while len(clades) > 1:
        new_clade = Clade(branch_length=random.random())
        children = random.sample(clades, 2)
        new_clade.clades.extend(children)
        clades = [c for c in clades if c not in children]
        clades.append(new_clade)
    tree.root = clades[0]
    return tree

def flat_tree(taxa):
    """Create a flat (star) tree with the given taxa."""
    tree = Tree(Clade())
    root_clade = Clade(branch_length=0)
    for taxon in taxa:
        clade = Clade(name=taxon, branch_length=random.random())
        root_clade.clades.append(clade)
    tree.root = root_clade
    return tree


# Hierarchical clustering as initial tree for MCMC sampling
def hierarchical_init_tree(cnv_data): 
    # Extract species names and CNV values
    species = list(cnv_data.keys())
    cnv_values = np.array(list(cnv_data.values()))

    # Compute pairwise distances
    distances = pdist(cnv_values, metric='euclidean')

    # Perform hierarchical clustering
    linkage_matrix = linkage(distances, method='average')

    # Convert scipy tree to Bio.Phylo tree
    def scipy_to_biopytree(node, parent=None):
        if node.is_leaf():
            return Clade(branch_length=0.1, name=species[node.id])
        else:
            clade = Clade(branch_length=0.1)
            clade.clades.append(scipy_to_biopytree(node.left, clade))
            clade.clades.append(scipy_to_biopytree(node.right, clade))
            return clade

    scipy_tree = to_tree(linkage_matrix)
    bio_tree = Phylo.BaseTree.Tree(root=scipy_to_biopytree(scipy_tree))

    return bio_tree
