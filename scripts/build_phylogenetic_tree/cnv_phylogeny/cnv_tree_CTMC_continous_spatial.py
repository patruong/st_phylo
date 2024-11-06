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
from concurrent.futures import ThreadPoolExecutor

import matplotlib.pyplot as plt
#from tree_generators import bifurcating_tree, random_tree, flat_tree, hierarchical_init_tree # initial trees
#from tree_modifiers import random_modify, ensemble_tree_modify

# Cache for pairwise distances
distance_cache = {}

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

#'''
def random_modify(tree):
    """Randomly modify the tree structure, including topology changes using NNI."""
    
    new_tree = copy.deepcopy(tree)
    
    # Get all internal clades (excluding root)
    internal_clades = [clade for clade in new_tree.find_clades(order='level') if clade != new_tree.root and not clade.is_terminal()]
    
    if not internal_clades:
        return new_tree  # Can't modify topology if no internal nodes
    
    # Select a random internal clade to perform NNI
    clade_to_modify = random.choice(internal_clades)
    
    # Get the parent of the clade
    parent = find_parent(new_tree.root, clade_to_modify)
    if parent is None or parent == new_tree.root:
        return new_tree  # Can't perform NNI on root
    
    # Get the siblings of the clade
    siblings = [clade for clade in parent.clades if clade != clade_to_modify]
    if not siblings:
        return new_tree  # No siblings to swap with
    
    sibling = random.choice(siblings)
    
    # Now swap subtrees between clade_to_modify and sibling
    # Swap a child of clade_to_modify with a child of sibling
    if clade_to_modify.is_terminal() or sibling.is_terminal():
        return new_tree  # Can't swap if any is terminal
    
    # Get child clades
    clade_children = clade_to_modify.clades
    sibling_children = sibling.clades
    
    if not clade_children or not sibling_children:
        return new_tree  # Can't swap if no children
    
    # Select random child from each to swap
    clade_child = random.choice(clade_children)
    sibling_child = random.choice(sibling_children)
    
    # Perform the swap
    clade_to_modify.clades.remove(clade_child)
    sibling.clades.remove(sibling_child)
    clade_to_modify.clades.append(sibling_child)
    sibling.clades.append(clade_child)
    
    # Modify branch lengths
    for clade in new_tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length *= random.uniform(0.8, 1.2)
    
    return new_tree
#'''
'''
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
'''
def find_parent(root, target_clade):
    """Find the parent of a given clade in the tree."""
    stack = [root]
    while stack:
        node = stack.pop()
        for child in node.clades:
            if child == target_clade:
                return node
            else:
                if child.clades:
                    stack.append(child)
    return None  # Parent not found

def hierarchical_init_tree(cnv_data, alpha=0.5):
    """
    Initialize a phylogenetic tree using a custom distance metric that combines
    spatial and genomic distances, with weighting factor `alpha`.
    
    Parameters:
    - cnv_data: Dictionary where each key is a species name and each value is a dictionary
                with "cnv_data" (genomic CNV data) and "coordinates" (spatial coordinates).
    - alpha: Weighting factor for spatial vs genomic distance (0 <= alpha <= 1).
    
    Returns:
    - bio_tree: Bio.Phylo tree based on the custom distance matrix.
    """
    # Extract species names, CNV values, and coordinates
    species = list(cnv_data.keys())
    cnv_values = np.array([data["cnv_data"] for data in cnv_data.values()])
    coordinates = np.array([data["coordinates"] for data in cnv_data.values()])

    # Compute pairwise genomic distances (D_genomic)
    genomic_distances = pdist(cnv_values, metric='euclidean')

    # Compute pairwise spatial distances (D_spatial)
    spatial_distances = pdist(coordinates, metric='euclidean')

    # Combine distances using the weighting factor alpha
    combined_distances = alpha * spatial_distances + (1 - alpha) * genomic_distances

    # Perform hierarchical clustering on the combined distances
    linkage_matrix = linkage(combined_distances, method='average')

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

def calculate_transition_prob_vectorized(rate_params, time, start_states, end_states, epsilon=1e-10):
    """Calculate transition probability for vectors of start and end states."""
    mu, sigma = rate_params
    start_states = np.maximum(start_states, epsilon)  # Avoid zero or negative start states
    end_states = np.maximum(end_states, epsilon)      # Avoid zero or negative end states

    m = np.log(start_states) + (mu - 0.5 * sigma**2) * time
    s = max(sigma * np.sqrt(time), epsilon)           # Ensure `s` is not zero
    scale = np.maximum(np.exp(m), epsilon)            # Ensure `scale` is positive

    return lognorm.pdf(end_states, s, scale=scale)

def tree_distance(tree, pair):
    """Compute and cache the distance between two sequences in the tree."""
    if pair in distance_cache:
        return distance_cache[pair]
    
    taxon1, taxon2 = pair
    clade1 = next((clade for clade in tree.get_terminals() if clade.name == taxon1), None)
    clade2 = next((clade for clade in tree.get_terminals() if clade.name == taxon2), None)
    
    if clade1 is None or clade2 is None:
        raise ValueError(f"One or both taxa not found in the tree: {taxon1}, {taxon2}")
    
    mrca = find_common_ancestor(tree, clade1, clade2)
    distance = tree.distance(clade1, mrca) + tree.distance(clade2, mrca)
    
    distance_cache[pair] = distance
    return distance

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

def cnv_likelihood(tree, cnv_data, rate_params, epsilon=1e-10, sample_fraction=0.1):
    """
    Compute likelihood of the given tree based on a subset of taxa pairs for faster computation.
    `cnv_data` should contain only CNV values, not spatial coordinates.
    """
    likelihood = 0.0
    # Extract only CNV data values for likelihood computation
    cnv_values_only = {key: data["cnv_data"] for key, data in cnv_data.items()}
    
    all_pairs = list(itertools.combinations(cnv_values_only.keys(), 2))
    sampled_pairs = random.sample(all_pairs, int(len(all_pairs) * sample_fraction))  # Sample subset
    
    def pair_likelihood(pair):
        cnv1, cnv2 = np.array(cnv_values_only[pair[0]]), np.array(cnv_values_only[pair[1]])
        distance = max(tree_distance(tree, pair), epsilon)
        transition_probs = np.maximum(calculate_transition_prob_vectorized(rate_params, distance, cnv1, cnv2), epsilon)
        return np.sum(np.log(transition_probs))  # Log-likelihood for this pair
    
    #with ThreadPoolExecutor() as executor:
    #    pair_likelihoods = list(executor.map(pair_likelihood, sampled_pairs))
    pair_likelihoods = [pair_likelihood(pair) for pair in sampled_pairs]

    return sum(pair_likelihoods) * (1 / sample_fraction)  # Scale for consistency

def cnv_spatial_likelihood(tree, cnv_data, rate_params, beta=0.5, epsilon=1e-10, sample_fraction=0.1):
    """
    Compute likelihood of the given tree based on both CNV and spatial proximity.
    `cnv_data` should contain both CNV values and spatial coordinates.
    
    Parameters:
    - tree: Phylogenetic tree.
    - cnv_data: Dictionary with "cnv_data" (genomic CNV data) and "coordinates" (spatial coordinates).
    - rate_params: Tuple with parameters (mu, sigma) for CNV model.
    - beta: Weighting factor for spatial vs CNV likelihood (0 <= beta <= 1).
    - epsilon: Small constant to avoid log(0).
    - sample_fraction: Fraction of pairs to sample for likelihood calculation (for efficiency).
    """
    likelihood = 0.0
    cnv_values_only = {key: data["cnv_data"] for key, data in cnv_data.items()}
    coordinates_only = {key: data["coordinates"] for key, data in cnv_data.items()}
    
    all_pairs = list(itertools.combinations(cnv_values_only.keys(), 2))
    sampled_pairs = random.sample(all_pairs, int(len(all_pairs) * sample_fraction))
    
    def pair_likelihood(pair):
        cnv1, cnv2 = np.array(cnv_values_only[pair[0]]), np.array(cnv_values_only[pair[1]])
        coord1, coord2 = np.array(coordinates_only[pair[0]]), np.array(coordinates_only[pair[1]])
        
        # Calculate genomic likelihood using CNV
        distance = max(tree_distance(tree, pair), epsilon)
        transition_probs = np.maximum(calculate_transition_prob_vectorized(rate_params, distance, cnv1, cnv2), epsilon)
        genomic_likelihood = np.sum(np.log(transition_probs))
        
        # Calculate spatial likelihood component
        spatial_distance = np.linalg.norm(coord1 - coord2)
        spatial_likelihood = -beta * spatial_distance  # Negative spatial distance increases likelihood for closer pairs
        
        # Combine both components (spatially closer spots are weighted more)
        return (1 - beta) * genomic_likelihood + spatial_likelihood

    pair_likelihoods = [pair_likelihood(pair) for pair in sampled_pairs]
    return sum(pair_likelihoods) * (1 / sample_fraction)  # Scale for consistency

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

def mcmc_tree_sampling_spatial(cnv_data, iterations=5000, beta=0.5):
    """
    Perform MCMC sampling to find the best phylogenetic tree based on both CNV data and spatial coordinates.
    
    Parameters:
    - cnv_data: Dictionary with "cnv_data" and "coordinates".
    - iterations: Number of MCMC iterations.
    - beta: Weighting factor for spatial influence.
    """
    tree = hierarchical_init_tree(cnv_data)
    rate_params = (0.01, 0.1)  # Starting values for rate parameters
    current_likelihood = cnv_spatial_likelihood(tree, cnv_data, rate_params, beta)
    best_tree = tree
    best_likelihood = current_likelihood
    
    likelihood_values = [current_likelihood]

    with tqdm(total=iterations, desc="MCMC Sampling") as pbar:
        for i in range(iterations):
            new_tree = random_modify(tree)
            new_likelihood = cnv_spatial_likelihood(new_tree, cnv_data, rate_params, beta)
            
            # Acceptance criterion
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

def convert_data_frame_to_input_dict(df):
    """Convert DataFrame to dictionary format for CNV data."""
    cnv_dict = df.apply(lambda row: row.tolist(), axis=1).to_dict()
    return cnv_dict

######## test on growmeatissue data

file_path = "/Users/patricktruong/git/st_phylo/data/synthetic/spatial_data/st-genome_profile.tsv"
genome_profile = pd.read_csv(file_path, sep='\t')
genome_profile = genome_profile.set_index(genome_profile.columns[0]).transpose()
genome_profile = genome_profile.loc[~(genome_profile.iloc[:, 1:] == 1.0).all(axis=1)]
genome_profile_log = np.log1p(genome_profile) + 1e-10  # Add epsilon to prevent log(0)


df = pd.read_csv("/Users/patricktruong/git/st_phylo/data/synthetic/spatial_data/st-meta.tsv", sep="\t")
annotations = pd.read_csv("/Users/patricktruong/git/st_phylo/data/synthetic/spatial_data/st-annotation.tsv", sep="\t", header=None, names=["Spot", "annotation"])

def reannotate(annotations, df):
    annotations["x"] = df.x
    annotations["y"] = df.y
    x_mid = (annotations['x'].max() + annotations['x'].min()) / 2
    y_mid = (annotations['y'].max() + annotations['y'].min()) / 2

    def update_annotation(row):
        if row['annotation'] == 'benign':
            return 'benign'
        elif row['annotation'] == 'aberrant':
            if (row['x'] > x_mid and row['y'] > y_mid):
                return 'A1'
            elif (row['x'] < x_mid and row['y'] < y_mid):
                return 'A2'
        return row['annotation']  # Keep original if none of the conditions match

    annotations['new_annotation'] = annotations.apply(update_annotation, axis=1)
    return annotations

def convert_data_frame_to_input_dict_with_coordinates(df, annotations):
    """
    Convert DataFrame to a dictionary format for CNV data with added spatial coordinates.
    Each entry will include the CNV values along with the (x, y) spatial coordinates.
    """
    # Create a dictionary with CNV data as lists
    cnv_dict = df.apply(lambda row: row.tolist(), axis=1).to_dict()

    # Add spatial coordinates from annotations
    for spot_id in cnv_dict.keys():
        spot_name = spot_id.split("{")[0]  # Extract the original spot name
        if spot_name in annotations.index:
            # Get coordinates as tuple (x, y)
            coordinates = (annotations.loc[spot_name, "x"], annotations.loc[spot_name, "y"])
            # Append coordinates to the CNV data for this spot
            cnv_dict[spot_id] = {"cnv_data": cnv_dict[spot_id], "coordinates": coordinates}
    
    return cnv_dict

annotations = reannotate(annotations, df)
annotations.new_annotation.unique()
annotations = reannotate(annotations, df)
annotations["annotation"] = annotations['new_annotation']
annotations["spot_id"] = annotations["Spot"] + "{" + annotations["annotation"] + "}"

spot_to_spot_id = dict(zip(annotations["Spot"], annotations["spot_id"]))

genome_profile_log.rename(index=spot_to_spot_id, inplace=True)
cnv_data = convert_data_frame_to_input_dict(genome_profile_log) # Code here just to illustrate how to get df to dict

#cnv_data['Spot_44{A2}']


annotations.set_index("Spot", inplace=True)  # Ensure annotations is indexed by 'Spot' for easy lookup
cnv_data_with_coords = convert_data_frame_to_input_dict_with_coordinates(genome_profile_log, annotations)
best_tree, best_likelihood, final_rate_params, likelihood_values = mcmc_tree_sampling(cnv_data_with_coords, iterations=1000) 
best_tree, best_likelihood, final_rate_params, likelihood_values = mcmc_tree_sampling_spatial(cnv_data_with_coords, iterations=1000, beta = 0.5) 
plot_likelihood_curve(likelihood_values)
plot_tree(best_tree)





########## Profiler
import cProfile
import pstats

cProfile.run('mcmc_tree_sampling(cnv_data, iterations=1)', 'profile_stats')
p = pstats.Stats('profile_stats')
p.sort_stats('cumulative').print_stats(10)  # Show the top 10 time-consuming functions