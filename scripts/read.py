import scanpy as sc
import pandas as pd
import os 
import matplotlib.pyplot as plt
import json
from PIL import Image
import numpy as np
Image.MAX_IMAGE_PIXELS = None  # Remove image size limit


import scanpy as sc
import pandas as pd

# Read the main data
adata = sc.read_visium('/home/ptruong/git/st_phylo/data/Count_matrices/Patient 1/Visium_with_annotation/H1_2/', 
                       count_file='filtered_feature_bc_matrix.h5',
                       load_images=True)

# Read the additional annotation file
annotations = pd.read_csv('path/to/H1_2_Final_Consensus_Annotations.csv')

# Add annotations to adata
adata.obs = adata.obs.join(annotations.set_index('barcode'), how='left')

# Print info about the loaded data
print(adata)

############################


def read_in_spatial(path, h5_file, tissue_position_file, scalefactors_file, image_file, sample_id, resize_img=True):
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

    # Read and add the image
    image_path = os.path.join(path, image_file)
    file_size = os.path.getsize(image_path)
    max_size = 1000 * 1024 * 1024  # 1000 MB, adjust as needed


    if file_size > max_size:
        print(f"Image file is too large ({file_size / (1024*1024):.2f} MB). Skipping.")
    else:
        if resize_img == True:
            with Image.open(image_path) as img:
                # Calculate new size (e.g., reduce to 25% of original size)
                new_size = tuple(dim // 4 for dim in img.size)
                img_resized = img.resize(new_size, Image.LANCZOS)
                img_array = np.array(img_resized)/255 #img_array_lowres
        else:
            img = Image.open(image_path)
            img_array = np.array(img)/255
        
    adata.uns['spatial'][sample_id]["images"] = {"hires": img_array}
    # add lowres
    # adata.uns['spatial'][sample_id]["images"] = {"lowres": img_array_lowres}

    print(adata.obsm['spatial'].shape)
    return adata


def plot_adata(adata):
    plt.scatter(adata.obs['array_col'], adata.obs['array_row'])
    plt.xlabel('Array Column')
    plt.ylabel('Array Row')
    plt.title('Spatial Distribution of Spots')
    plt.show()

if __name__ == "__main__":
    path = "/home/ptruong/git/st_phylo/data/Count_matrices/Patient 1/Visium_with_annotation/H1_2/"
    h5_file = "filtered_feature_bc_matrix.h5"
    tissue_position_file = "tissue_positions_list.csv"
    scalefactors_file = "scalefactors_json.json"
    image_file = "H1_2.jpg"
    sample_id = "Patient_2_H1_2"

    adata = read_in_spatial(path, h5_file, tissue_position_file, scalefactors_file, image_file, sample_id, resize_img=True)

    plot_adata(adata)