import numpy as np
import itertools
from Bio import Phylo
from io import StringIO
import random
import copy


# Simple Jukes-Cantor model likelihood function
def jukes_cantor_likelihood(tree, sequences, epsilon=1e-10):
    """Compute likelihood of the given tree under the Jukes-Cantor model."""
    likelihood = 0.0
    
    # Iterate over all pairs of sequences
    for pair in itertools.combinations(sequences.keys(), 2):
        seq1, seq2 = sequences[pair[0]], sequences[pair[1]]
        distance = tree_distance(tree, pair)
        
        # Iterate over each site in the alignment
        site_likelihood = 0.0
        for a, b in zip(seq1, seq2):
            if a == b:  # No difference
                p_same = 0.25 + 0.75 * np.exp(-4/3 * distance)
                site_likelihood += np.log(max(p_same, epsilon))  # Avoid log(0)
            else:       # Difference
                p_diff = 0.25 * (1 - np.exp(-4/3 * distance))
                site_likelihood += np.log(max(p_diff, epsilon))  # Avoid log(0)
        
        # Add to total likelihood
        likelihood += site_likelihood

    return likelihood

# Define a simple function to simulate tree distances
def tree_distance(tree, pair):
    """Compute the distance between two sequences in the tree."""
    taxon1, taxon2 = pair[0], pair[1]
    
    # Find the clades corresponding to the taxa names
    clade1 = tree.find_any(name=taxon1)
    clade2 = tree.find_any(name=taxon2)
    
    # Calculate the distance between the two clades
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

def random_modify(tree):
    """Randomly modify the tree structure."""
    new_tree = copy.deepcopy(tree)
    if len(list(new_tree.find_clades())) <= 3:
        return new_tree
    
    non_root_clades = list(new_tree.find_clades(terminal=False))[1:]  # Exclude root
    if not non_root_clades:
        return new_tree
    
    clade_to_modify = random.choice(non_root_clades)
    if len(clade_to_modify.clades) >= 2:
        # Randomly rearrange the children of this clade
        random.shuffle(clade_to_modify.clades)
    
    # Randomly adjust branch lengths
    for clade in new_tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length *= random.uniform(0.8, 1.2)
    
    return new_tree

def mcmc_tree_sampling(sequences, iterations=10000):
    """Perform MCMC sampling to find the best phylogenetic tree."""
    taxa = list(sequences.keys())
    tree = random_tree(taxa)
    
    current_likelihood = jukes_cantor_likelihood(tree, sequences)
    best_tree = tree
    best_likelihood = current_likelihood

    for i in range(iterations):
        new_tree = random_modify(tree)
        new_likelihood = jukes_cantor_likelihood(new_tree, sequences)
        
        if new_likelihood > current_likelihood or random.random() < np.exp(new_likelihood - current_likelihood):
            tree = new_tree
            current_likelihood = new_likelihood

            if new_likelihood > best_likelihood:
                best_tree = new_tree
                best_likelihood = new_likelihood
    
    return best_tree, best_likelihood

sequences = {
    # Clade 1: Similar sequences
    "A": "ATCGTACGATCG",  # Species A
    "B": "ATCGTACGATGG",  # Species B
    "C": "ATCGTACGCTGG",  # Species C
    
    # Clade 2: Similar sequences but different from Clade 1
    "D": "TTGCTAGGGTCA",  # Species D
    "E": "TTGCTAGGGTGA",  # Species E
    "F": "TTGCTAGGATGA",  # Species F
    
    # Intermediate taxa to make the tree more complex
    "G": "ATCGTATGCTGG",  # Closer to Clade 1
    "H": "TTGCTACGGTGA",  # Closer to Clade 2
}

# Run MCMC to infer the best tree
best_tree, best_likelihood = mcmc_tree_sampling(sequences, iterations=10000)

# Print the best tree
tree_str = StringIO()
Phylo.write(best_tree, tree_str, "newick")
print(f"Best tree (Newick format): {tree_str.getvalue()}")
print(f"Best likelihood: {best_likelihood}")

# Visualize the best tree
Phylo.draw(best_tree)
