import scanpy as sc
import pandas as pd
import os 
import matplotlib.pyplot as plt
import json

def read_in_spatial(path, h5_file, tissue_position_file, scalefactors_file, sample_id):
    file_path = os.path.join(path, h5_file)
    adata = sc.read_10x_h5(file_path)
    adata.var_names_make_unique()

    # Read positions file
    positions_path = os.path.join(path, tissue_position_file)
    positions = pd.read_csv(positions_path, header=None, names=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres'])

    # Set barcode as index
    positions.set_index('barcode', inplace=True)

    # Filter positions to match adata
    positions_filtered = positions.loc[adata.obs_names]

    # Now assign to adata.obsm
    adata.obsm['spatial'] = positions_filtered[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values

    # You can also add these as columns in adata.obs if you want
    adata.obs['array_row'] = positions_filtered['array_row']
    adata.obs['array_col'] = positions_filtered['array_col']
    adata.obs['in_tissue'] = positions_filtered['in_tissue']

    # Read and add scalefactors
    scalefactors_path = os.path.join(path, scalefactors_file)
    with open(scalefactors_path, 'r') as f:
        scalefactors = json.load(f)
    
    # Initialize the 'spatial' key in adata.uns if it doesn't exist
    if 'spatial' not in adata.uns:
        adata.uns['spatial'] = {}
    
    # Add the scalefactors to adata.uns
    adata.uns['spatial'][sample_id] = {"scalefactors": scalefactors}

    print(adata.obsm['spatial'].shape)
    return adata


def plot_adata(adata):
    plt.scatter(adata.obs['array_col'], adata.obs['array_row'])
    plt.xlabel('Array Column')
    plt.ylabel('Array Row')
    plt.title('Spatial Distribution of Spots')
    plt.show()

# Usage
path = "/home/ptruong/git/st_phylo/data/Count_matrices/Patient 1/Visium_with_annotation/H1_2/"
h5_file = "filtered_feature_bc_matrix.h5"
tissue_position_file = "tissue_positions_list.csv"
scalefactors_file = "scalefactors_json.json"
sample_id = "Patient_2_H2_2"

adata2 = read_in_spatial(path, h5_file, tissue_position_file, scalefactors_file, sample_id)

plot_adata(adata)