import pandas as pd
import numpy as np
import itertools
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

from scipy.stats import lognorm
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import pdist
from scipy.optimize import minimize

import random
import copy
from tqdm import tqdm

import matplotlib.pyplot as plt
from tree_generators import bifurcating_tree, random_tree, flat_tree, hierarchical_init_tree # initial trees
from tree_modifiers import random_modify, ensemble_tree_modify


def plot_tree(tree):
    """Plot the phylogenetic tree using matplotlib."""
    fig, ax = plt.subplots(figsize=(10, 8))
    Phylo.draw(tree, axes=ax, do_show=False)
    
    # Customize the plot
    ax.set_xlabel('Branch Length')
    ax.set_ylabel('Taxa')
    ax.set_title('Phylogenetic Tree')
    
    # Remove the default axis on the right and top of the plot
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Show the plot
    plt.tight_layout()
    plt.show()

def plot_likelihood_curve(likelihood_values):
    """Plot the likelihood curve from MCMC sampling."""
    plt.figure(figsize=(12, 6))
    plt.plot(likelihood_values)
    plt.title('MCMC Sampling Likelihood Curve')
    plt.xlabel('Iteration')
    plt.ylabel('Log-Likelihood')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()


def cnv_likelihood(tree, cnv_data, rate_params, epsilon=1e-10):
    """Compute likelihood of the given tree based on continuous CNV data using log-normal distribution."""
    likelihood = 0.0
    for pair in itertools.combinations(cnv_data.keys(), 2):
        cnv1, cnv2 = cnv_data[pair[0]], cnv_data[pair[1]]
        distance = tree_distance(tree, pair)
        
        for a, b in zip(cnv1, cnv2):
            transition_prob = calculate_transition_prob(rate_params, distance, a, b)
            likelihood += np.log(max(transition_prob, epsilon))
    
    return likelihood

def calculate_transition_prob(rate_params, time, start_state, end_state):
    """Calculate transition probability using log-normal distribution."""
    mu, sigma = rate_params
    # Adjust parameters for log-normal distribution
    m = np.log(start_state) + (mu - 0.5 * sigma**2) * time
    s = sigma * np.sqrt(time)
    return lognorm.pdf(end_state, s, scale=np.exp(m))

#def tree_distance(tree, pair):
#    """Compute the distance between two sequences in the tree."""
#    taxon1, taxon2 = pair[0], pair[1]
#    clade1 = tree.find_any(name=taxon1)
#    clade2 = tree.find_any(name=taxon2)
 #   dist = tree.distance(clade1, clade2)
 #   return dist

def find_common_ancestor(tree, taxon1, taxon2):
    """Find the most recent common ancestor of two taxa in the tree."""
    path1 = tree.get_path(taxon1)
    path2 = tree.get_path(taxon2)
    
    # Find the last common element in both paths
    for i, (node1, node2) in enumerate(zip(path1, path2)):
        if node1 != node2:
            return path1[i-1]
    
    # If one path is a prefix of the other, return the last element of the shorter path
    return path1[-1] if len(path1) < len(path2) else path2[-1]

def tree_distance(tree, pair):
    """Compute the distance between two sequences in the tree."""
    taxon1, taxon2 = pair[0], pair[1]
    
    # Find the corresponding Clade objects for the taxa
    clade1 = next((clade for clade in tree.get_terminals() if clade.name == taxon1), None)
    clade2 = next((clade for clade in tree.get_terminals() if clade.name == taxon2), None)
    
    if clade1 is None or clade2 is None:
        raise ValueError(f"One or both taxa not found in the tree: {taxon1}, {taxon2}")
    
    # Find the most recent common ancestor
    mrca = find_common_ancestor(tree, clade1, clade2)
    
    # Calculate the distance as the sum of branch lengths from each taxon to the MRCA
    distance = (tree.distance(clade1, mrca) + tree.distance(clade2, mrca))
    
    return distance


def mcmc_tree_sampling(cnv_data, iterations=5000):
    """Perform MCMC sampling to find the best phylogenetic tree based on continuous CNV data."""
    tree = hierarchical_init_tree(cnv_data)
    rate_params = (0.01, 0.1)  # Example values, should be tuned
    current_likelihood = cnv_likelihood(tree, cnv_data, rate_params)
    best_tree = tree
    best_likelihood = current_likelihood

    # List to store likelihood values for plotting
    likelihood_values = [current_likelihood]

    with tqdm(total=iterations, desc="MCMC Sampling") as pbar:
        for i in range(iterations):
            new_tree = random_modify(tree)
            new_likelihood = cnv_likelihood(new_tree, cnv_data, rate_params)
            if new_likelihood > current_likelihood or random.random() < np.exp(new_likelihood - current_likelihood):
                tree = new_tree
                current_likelihood = new_likelihood
                if new_likelihood > best_likelihood:
                    best_tree = new_tree
                    best_likelihood = new_likelihood
            
            likelihood_values.append(current_likelihood)
            
            if i % 100 == 0:
                rate_params = update_rate_params(rate_params, tree, cnv_data)
            
            pbar.update(1)
            if i % 100 == 0:
                pbar.set_description(f"MCMC Sampling (Best Likelihood: {best_likelihood:.2f})")

    return best_tree, best_likelihood, rate_params, likelihood_values


# Random update rate for testing
#def update_rate_params(current_params, tree, cnv_data):
#    """Update rate parameters based on current tree and data."""
#    # This is a placeholder. In practice, you'd implement a method to optimize these parameters
#    # based on the current tree structure and observed data.
#    mu, sigma = current_params
#    mu *= random.uniform(0.9, 1.1)
#    sigma *= random.uniform(0.9, 1.1)
#    return (mu, sigma)


def update_rate_params(current_params, tree, cnv_data):
    """Optimize rate parameters (mu, sigma) to maximize likelihood based on current tree and CNV data."""
    
    # Define the objective function as the negative likelihood
    def objective(params):
        mu, sigma = params
        if sigma <= 0:  # Ensure sigma is positive
            return np.inf  # Penalize non-positive sigma
        return -cnv_likelihood(tree, cnv_data, (mu, sigma))  # Negative for minimization
    
    # Set bounds for mu and sigma to prevent unrealistic values
    bounds = [(1e-6, 5.0), (1e-6, 5.0)]  # Adjust bounds based on expected parameter ranges
    
    # Run the optimizer starting from the current parameters
    result = minimize(objective, current_params, bounds=bounds, method='L-BFGS-B')
    
    # Return the optimized parameters if successful, otherwise return the current ones
    if result.success:
        return result.x  # Optimized mu, sigma
    else:
        print("Optimization failed; returning current parameters.")
        return current_params  # Fallback in case optimization fails

cnv_data = {
    # Clade A
    'species1_A': [10.3, 1.0, 1.0, 1.0, 1.0],
    'species2_A': [10.2, 1.0, 1.0, 1.0, 1.0],
    'species3_A': [10.1, 1.0, 1.0, 1.0, 1.0],
    'species4_A': [10.4, 1.0, 1.0, 1.0, 1.0],
    'species5_A': [10.0, 1.0, 1.0, 1.0, 1.0],
    'species6_A': [10.3, 1.0, 1.0, 1.0, 1.0],

    # Clade B
    'species7_B': [1.0, 10.5, 1.0, 1.0, 1.0],
    'species8_B': [1.0, 10.4, 1.0, 1.0, 1.0],
    'species9_B': [1.0, 10.2, 1.0, 1.0, 1.0],
    'species10_B': [1.0, 10.3, 1.0, 1.0, 1.0],
    'species11_B': [1.0, 10.6, 1.0, 1.0, 1.0],
    'species12_B': [1.0, 10.1, 1.0, 1.0, 1.0],

    # Clade C
    'species13_C': [1.0, 1.0, 1.0, 9.0, 10.0],
    'species14_C': [1.0, 1.0, 1.0, 9.2, 10.1],
    'species15_C': [1.0, 1.0, 1.0, 9.1, 10.2],
    'species16_C': [1.0, 1.0, 1.0, 8.9, 9.9]
}


# Run the MCMC sampling
best_tree, best_likelihood, final_rate_params, likelihood_values = mcmc_tree_sampling(cnv_data, iterations=1000)

print(f"Best likelihood: {best_likelihood}")
print(f"Final rate parameters: mu = {final_rate_params[0]}, sigma = {final_rate_params[1]}")

# Plot
plot_likelihood_curve(likelihood_values)
plot_tree(best_tree)

# Make this work for this data
cnv_data = pd.DataFrame(cnv_data) #for vis
cnv_data = cnv_data.T
np.array(cnv_data)


