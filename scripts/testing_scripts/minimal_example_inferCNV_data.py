import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import seaborn as sns
import matplotlib.pyplot as plt
from PIL import Image
Image.MAX_IMAGE_PIXELS = None  # Remove image size limit
import os
os.chdir("/home/ptruong/git/st_phylo/test")

img_file = "H3_2.jpg"
image_path = os.path.join(path, img_file)

#with Image.open(image_path) as img:
#    img_array = np.array(img)
    
# Open the image
with Image.open(image_path) as img:
    constant = 2.02
    scaling_factor_hires = 0.04111842 * constant # from scaling_factors_json.json
    scaling_factor_lowres = 0.01233552 * constant # from scaling_factors_json.json

    # hires = 0.04111842
    # lowres = 0.01233552
    new_width_hires = int(img.width * scaling_factor_hires)
    new_height_hires = int(img.height * scaling_factor_hires)
    img_hires = img.resize((new_width_hires, new_height_hires), Image.LANCZOS)
    new_width_lowres = int(img.width * scaling_factor_lowres)
    new_height_lowres = int(img.height * scaling_factor_lowres)
    img_lowres = img.resize((new_width_lowres, new_height_lowres), Image.LANCZOS)

    img_hires.save(os.path.join(path, "spatial/tissue_hires_image.png"), format='PNG')
    img_lowres.save(os.path.join(path, "spatial/tissue_lowres_image.png"), format='PNG')


#adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
path = "/home/ptruong/git/st_phylo/test"

#adata = sc.read_visium(path, load_images = False)
adata = sc.read_visium(path)
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
sc.pl.spatial(adata, img_key="lowres", color=["total_counts", "n_genes_by_counts"])