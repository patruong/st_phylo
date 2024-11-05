from cnv_tree import *

# Set a random seed for reproducibility
np.random.seed(42)

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
best_tree, best_likelihood = mcmc_tree_sampling(cnv_data_pure, iterations=10000)

# Print the best tree
tree_str = StringIO()
Phylo.write(best_tree, tree_str, "newick")
print(f"Best tree (Newick format): {tree_str.getvalue()}")
print(f"Best likelihood: {best_likelihood}")

# Visualize the best tree
Phylo.draw(best_tree)


