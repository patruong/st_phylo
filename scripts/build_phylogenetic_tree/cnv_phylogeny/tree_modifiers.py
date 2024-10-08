import random
import copy
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

def ensemble_tree_modify(tree, cnv_data):
    """
    Modify the tree structure using biologically relevant operations while preserving all taxa.
    
    This improved tree modification method offers several advantages:
    1. Biologically relevant operations: It implements three common phylogenetic tree modification operations:
       - Nearest Neighbor Interchange (NNI)
       - Subtree Pruning and Regrafting (SPR)
       - Branch length adjustment
    2. Weighted operation selection: The method uses weighted probabilities to choose between different operations.
    3. Safe modifications: Ensures that all taxa are preserved during modifications.
    """
    new_tree = copy.deepcopy(tree)
    
    # List of possible tree modification operations
    operations = [
        (0.35, safe_nearest_neighbor_interchange),
        (0.35, safe_subtree_pruning_and_regrafting),
        (0.3, adjust_branch_lengths),
    ]
    
    # Choose an operation based on weighted probabilities
    operation = random.choices(operations, weights=[op[0] for op in operations])[0][1]
    
    # Apply the chosen operation
    modified_tree = operation(new_tree)
    
    # Ensure all taxa are present
    all_taxa = set(cnv_data.keys())
    tree_taxa = set(leaf.name for leaf in modified_tree.get_terminals())
    
    if all_taxa != tree_taxa:
        print(f"Warning: Taxa mismatch. Reverting to original tree.")
        return tree
    
    return modified_tree

def safe_nearest_neighbor_interchange(tree):
    """Perform Nearest Neighbor Interchange (NNI) on the tree, ensuring all taxa are preserved."""
    internal_nodes = [clade for clade in tree.find_clades(terminal=False) if len(clade.clades) > 2]
    if not internal_nodes:
        return tree
    
    node = random.choice(internal_nodes)
    if len(node.clades) < 3:
        return tree
    
    i, j = random.sample(range(len(node.clades)), 2)
    node.clades[i], node.clades[j] = node.clades[j], node.clades[i]
    
    return tree

def find_parent(tree, node):
    """Find the parent of a given node in the tree."""
    for clade in tree.find_clades():
        if node in clade.clades:
            return clade
    return None

def safe_subtree_pruning_and_regrafting(tree):
    """Perform Subtree Pruning and Regrafting (SPR) on the tree, ensuring all taxa are preserved."""
    internal_nodes = list(tree.find_clades(terminal=False))
    if len(internal_nodes) < 2:
        return tree
    
    # Choose a random internal node to prune (excluding the root)
    prune_node = random.choice(internal_nodes[1:])
    parent = find_parent(tree, prune_node)
    
    if parent is None or parent == tree.root:
        return tree
    
    # Remove the pruned node from its parent
    parent.clades.remove(prune_node)
    
    # Choose a random place to regraft (excluding the parent and the pruned node itself)
    possible_regraft_nodes = [node for node in internal_nodes if node != prune_node and node != parent]
    
    if not possible_regraft_nodes:
        # If no suitable regraft location, reattach to the original parent
        parent.clades.append(prune_node)
        return tree
    
    regraft_node = random.choice(possible_regraft_nodes)
    regraft_node.clades.append(prune_node)
    
    return tree

def adjust_branch_lengths(tree):
    """Adjust branch lengths of the tree."""
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length *= random.uniform(0.8, 1.2)
    return tree

# The following function is kept for reference but not used in the safe version
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