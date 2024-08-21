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
    scaling_factor_hires = 0.04111842 * constant
    scaling_factor_lowres = 0.01233552 * constant

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


fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1],
)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3],
)

plt.show()

sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(
    adata, key_added="clusters", directed=False, n_iterations=2
)



#plt.rcParams["figure.figsize"] = (4, 4)
#sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
sc.pl.spatial(adata, img_key="lowres", color=["total_counts", "n_genes_by_counts"])


np.shape(adata.uns["spatial"]['V10L27-044_B1']["images"]["hires"])

np.shape(adata.uns["spatial"]['V10L27-044_B1']["images"]["lowres"])

sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)

sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    groups=["5", "9"],
    crop_coord=[7000, 10000, 0, 6000],
    alpha=0.5,
    size=1.3,
)

sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="9", n_genes=10, groupby="clusters")

sc.pl.spatial(adata, img_key="hires", color=["clusters", "CR2"])

sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)

