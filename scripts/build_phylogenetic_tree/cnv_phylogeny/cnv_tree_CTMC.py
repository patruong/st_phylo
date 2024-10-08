import numpy as np
import itertools
from Bio import Phylo
from scipy.linalg import expm
import random
import copy


"""
In this updated model:

States represent average copy number (CN) ranges rather than discrete states.
Transitions represent gradual changes in the average CN within a spot.
λ_inc represents the rate of increasing average CN.
μ_dec represents the rate of decreasing average CN.

"""

def cnv_likelihood(tree, cnv_data, rate_matrix, epsilon=1e-10):
    """Compute likelihood of the given tree based on CNV data using a continuous-time Markov chain model."""
    likelihood = 0.0
    for pair in itertools.combinations(cnv_data.keys(), 2):
        cnv1, cnv2 = cnv_data[pair[0]], cnv_data[pair[1]]
        distance = tree_distance(tree, pair)
        
        for a, b in zip(cnv1, cnv2):
            transition_prob = calculate_transition_prob(rate_matrix, distance, a, b)
            likelihood += np.log(max(transition_prob, epsilon))
    
    return likelihood

def calculate_transition_prob(rate_matrix, time, start_state, end_state):
    """Calculate transition probability using matrix exponentiation."""
    transition_matrix = expm(rate_matrix * time)
    return transition_matrix[start_state, end_state]

def initialize_rate_matrix(n_states):
    """Initialize a simple rate matrix for the continuous-time Markov chain."""
    rate_matrix = np.zeros((n_states, n_states))
    for i in range(n_states):
        if i > 0:
            rate_matrix[i, i-1] = 1.0  # Rate of deletion
        if i < n_states - 1:
            rate_matrix[i, i+1] = 1.0  # Rate of duplication
        rate_matrix[i, i] = -np.sum(rate_matrix[i, :])
    return rate_matrix

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

def random_modify(tree):
    """Randomly modify the tree structure."""
    new_tree = copy.deepcopy(tree)
    if len(list(new_tree.find_clades())) <= 3:
        return new_tree
    non_root_clades = list(new_tree.find_clades(terminal=False))[1:]
    if not non_root_clades:
        return new_tree
    clade_to_modify = random.choice(non_root_clades)
    if len(clade_to_modify.clades) >= 2:
        random.shuffle(clade_to_modify.clades)
    for clade in new_tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length *= random.uniform(0.8, 1.2)
    return new_tree

def mcmc_tree_sampling(cnv_data, iterations=10000):
    """Perform MCMC sampling to find the best phylogenetic tree based on CNV data."""
    taxa = list(cnv_data.keys())
    tree = random_tree(taxa)
    
    # Initialize the rate matrix
    max_cnv = max(max(cnv) for cnv in cnv_data.values())
    rate_matrix = initialize_rate_matrix(max_cnv + 1)
    
    current_likelihood = cnv_likelihood(tree, cnv_data, rate_matrix)
    best_tree = tree
    best_likelihood = current_likelihood
    
    for i in range(iterations):
        new_tree = random_modify(tree)
        new_likelihood = cnv_likelihood(new_tree, cnv_data, rate_matrix)
        
        if new_likelihood > current_likelihood or random.random() < np.exp(new_likelihood - current_likelihood):
            tree = new_tree
            current_likelihood = new_likelihood
            
        if new_likelihood > best_likelihood:
            best_tree = new_tree
            best_likelihood = new_likelihood
    
    return best_tree, best_likelihood

# Example usage
cnv_data = {
    'species1': [2, 3, 2, 4],
    'species2': [2, 2, 3, 3],
    'species3': [3, 3, 2, 2],
    'species4': [2, 4, 3, 2]
}

max_cnv = max(max(cnv) for cnv in cnv_data.values())
rate_matrix = initialize_rate_matrix(max_cnv + 1)

best_tree, best_likelihood = mcmc_tree_sampling(cnv_data)
print(f"Best likelihood: {best_likelihood}")
Phylo.draw_ascii(best_tree)

