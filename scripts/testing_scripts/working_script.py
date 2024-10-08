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


### ToDo we need to format all the data in a way so that i can be parsed with sc.read_visium

