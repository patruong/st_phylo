import numpy as np
import itertools
from Bio import Phylo
from io import StringIO
import random
import copy


def cnv_likelihood(tree, cnv_data, epsilon=1e-10):
    """Compute likelihood of the given tree based on CNV data."""
    likelihood = 0.0
    
    # Iterate over all pairs of species
    for pair in itertools.combinations(cnv_data.keys(), 2):
        cnv1, cnv2 = cnv_data[pair[0]], cnv_data[pair[1]]
        distance = tree_distance(tree, pair)
        
        # Compute likelihood based on CNV differences
        for a, b in zip(cnv1, cnv2):
            diff = abs(a - b)
            # Placeholder: simple Gaussian-like likelihood
            site_likelihood = np.exp(-diff**2 / (2 * distance))
            likelihood += np.log(max(site_likelihood, epsilon))  # Avoid log(0)

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

# Update the mcmc_tree_sampling function to use cnv_likelihood
def mcmc_tree_sampling(cnv_data, iterations=10000):
    """Perform MCMC sampling to find the best phylogenetic tree based on CNV data."""
    taxa = list(cnv_data.keys())
    tree = random_tree(taxa)
    
    current_likelihood = cnv_likelihood(tree, cnv_data)
    best_tree = tree
    best_likelihood = current_likelihood

    for i in range(iterations):
        new_tree = random_modify(tree)
        new_likelihood = cnv_likelihood(new_tree, cnv_data)
        
        if new_likelihood > current_likelihood or random.random() < np.exp(new_likelihood - current_likelihood):
            tree = new_tree
            current_likelihood = new_likelihood

            if new_likelihood > best_likelihood:
                best_tree = new_tree
                best_likelihood = new_likelihood
    
    return best_tree, best_likelihood

cnv_data_pure = {
    # Clade 1: Similar CNV profiles
    "A": [1.0, 1.5, 1.0, 0.5, 1.0, 1.0, 2.0],  # Species A
    "B": [1.0, 1.5, 1.0, 0.5, 1.0, 1.0, 1.5],  # Species B
    "C": [1.0, 1.5, 1.0, 0.5, 1.0, 0.5, 1.5],  # Species C
    
    # Clade 2: Similar CNV profiles but different from Clade 1
    "D": [1.5, 1.5, 0.5, 1.0, 1.5, 2.0, 1.0],  # Species D
    "E": [1.5, 1.5, 0.5, 1.0, 1.5, 2.0, 0.5],  # Species E
    "F": [1.5, 1.5, 0.5, 1.0, 1.5, 0.5, 0.5],  # Species F
    
    # Intermediate taxa to make the tree more complex
    "G": [1.0, 1.5, 1.0, 0.5, 1.0, 1.5, 1.5],  # Closer to Clade 1
    "H": [1.5, 1.5, 0.5, 1.0, 1.0, 2.0, 0.5],  # Closer to Clade 2
}

# Set a random seed for reproducibility
np.random.seed(42)

def add_noise(value, noise_level=0.1):
    """Add Gaussian noise to a value."""
    return value + np.random.normal(0, noise_level * value)

cnv_data = {
    # Clade 1: Similar CNV profiles
    "A": [add_noise(1.0), add_noise(1.5), add_noise(1.0), add_noise(0.5), add_noise(1.0), add_noise(1.0), add_noise(2.0)],  # Species A
    "B": [add_noise(1.0), add_noise(1.5), add_noise(1.0), add_noise(0.5), add_noise(1.0), add_noise(1.0), add_noise(1.5)],  # Species B
    "C": [add_noise(1.0), add_noise(1.5), add_noise(1.0), add_noise(0.5), add_noise(1.0), add_noise(0.5), add_noise(1.5)],  # Species C
    
    # Clade 2: Similar CNV profiles but different from Clade 1
    "D": [add_noise(1.5), add_noise(1.5), add_noise(0.5), add_noise(1.0), add_noise(1.5), add_noise(2.0), add_noise(1.0)],  # Species D
    "E": [add_noise(1.5), add_noise(1.5), add_noise(0.5), add_noise(1.0), add_noise(1.5), add_noise(2.0), add_noise(0.5)],  # Species E
    "F": [add_noise(1.5), add_noise(1.5), add_noise(0.5), add_noise(1.0), add_noise(1.5), add_noise(0.5), add_noise(0.5)],  # Species F
    
    # Intermediate taxa to make the tree more complex
    "G": [add_noise(1.0), add_noise(1.5), add_noise(1.0), add_noise(0.5), add_noise(1.0), add_noise(1.5), add_noise(1.5)],  # Closer to Clade 1
    "H": [add_noise(1.5), add_noise(1.5), add_noise(0.5), add_noise(1.0), add_noise(1.0), add_noise(2.0), add_noise(0.5)],  # Closer to Clade 2
}

# Print the continuous CNV data
for species, cnv_values in cnv_data.items():
    print(f"{species}: {[round(value, 3) for value in cnv_values]}")


##### Here this wont work because it is not yet adapter for continuous data.
# Run MCMC to infer the best tree
best_tree, best_likelihood = mcmc_tree_sampling(cnv_data, iterations=10000)

# Print the best tree
tree_str = StringIO()
Phylo.write(best_tree, tree_str, "newick")
print(f"Best tree (Newick format): {tree_str.getvalue()}")
print(f"Best likelihood: {best_likelihood}")

# Visualize the best tree
Phylo.draw(best_tree)


