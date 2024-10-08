import numpy as np
import pymc as pm
import arviz as az
from Bio import Phylo
from io import StringIO
import random
import itertools

# Helper functions (tree_distance, random_tree, random_modify) remain the same
def tree_distance(tree, pair):
    """Compute the distance between two sequences in the tree."""
    taxon1, taxon2 = pair[0], pair[1]
    clade1 = tree.find_any(name=taxon1)
    clade2 = tree.find_any(name=taxon2)
    dist = tree.distance(clade1, clade2)
    return dist

def random_tree(taxa):
    """Create a random binary tree with the given taxa."""
    tree = Phylo.BaseTree.Tree(Phylo.BaseTree.Clade())
    clades = [Phylo.BaseTree.Clade(name=taxon, branch_length=random.random()) for taxon in taxa]
    while len(clades) > 1:
        new_clade = Phylo.BaseTree.Clade(branch_length=random.random())
        children = random.sample(clades, 2)
        new_clade.clades.extend(children)
        clades = [c for c in clades if c not in children]
        clades.append(new_clade)
    tree.root = clades[0]
    return tree

def jukes_cantor_likelihood(tree, sequences):
    """Compute log-likelihood of the given tree under the Jukes-Cantor model."""
    log_likelihood = 0.0
    
    for pair in itertools.combinations(sequences.keys(), 2):
        seq1, seq2 = sequences[pair[0]], sequences[pair[1]]
        distance = tree_distance(tree, pair)
        
        for a, b in zip(seq1, seq2):
            if a == b:
                p_same = 0.25 + 0.75 * np.exp(-4/3 * distance)
                log_likelihood += np.log(p_same)
            else:
                p_diff = 0.25 * (1 - np.exp(-4/3 * distance))
                log_likelihood += np.log(p_diff)
    
    return log_likelihood

def pymc_tree_inference(sequences, n_samples=1000, tune=1000):
    """Perform Bayesian phylogenetic inference using PyMC."""
    taxa = list(sequences.keys())
    n_taxa = len(taxa)
    
    with pm.Model() as model:
        # Prior on tree topology (uniform over all possible trees)
        tree = pm.Categorical("tree", p=np.ones(n_taxa-2) / (n_taxa-2))
        
        # Prior on branch lengths (exponential distribution)
        branch_lengths = pm.Exponential("branch_lengths", lam=10, shape=2*n_taxa-3)
        
        # Likelihood function
        def likelihood_func(tree, branch_lengths):
            # Convert discrete tree representation and branch lengths to a Phylo tree
            phylo_tree = random_tree(taxa)  # Start with a random tree structure
            for i, clade in enumerate(phylo_tree.find_clades()):
                if i < len(branch_lengths):
                    clade.branch_length = branch_lengths[i]
            
            return jukes_cantor_likelihood(phylo_tree, sequences)
        
        likelihood = pm.Deterministic("likelihood", likelihood_func(tree, branch_lengths))
        
        # Sample from the posterior
        trace = pm.sample(n_samples, tune=tune, chains=4, cores=4)
    
    return trace, model

# Example usage
sequences = {
    "A": "ATCGTACGATCG",
    "B": "ATCGTACGATGG",
    "C": "ATCGTACGCTGG",
    "D": "TTGCTAGGGTCA",
    "E": "TTGCTAGGGTGA",
    "F": "TTGCTAGGATGA",
    "G": "ATCGTATGCTGG",
    "H": "TTGCTACGGTGA",
}

# Run Bayesian inference
trace, model = pymc_tree_inference(sequences)

# Analyze results
summary = az.summary(trace)
print(summary)

# Extract the maximum a posteriori (MAP) estimate
map_estimate = pm.find_MAP(model=model)

# Convert MAP estimate to a Phylo tree
map_tree = random_tree(list(sequences.keys()))
for i, clade in enumerate(map_tree.find_clades()):
    if i < len(map_estimate['branch_lengths']):
        clade.branch_length = map_estimate['branch_lengths'][i]

# Print the best tree
tree_str = StringIO()
Phylo.write(map_tree, tree_str, "newick")
print(f"Best tree (Newick format): {tree_str.getvalue()}")

# Visualize the best tree
Phylo.draw(map_tree)