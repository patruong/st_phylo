import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
import scanpy as sc
import numpy as np
from anndata import AnnData
import squidpy as sq
from PIL import Image


def plot_genome_profile_heatmap(file_path, clustering_method='ward', cmap="viridis"):
    """
    Plots a genome profile heatmap with hierarchical clustering on spots.

    Parameters:
    - file_path (str): Path to the genome profile TSV file.
    - clustering_method (str): Clustering method to use for hierarchical clustering (default is 'ward').
    - cmap (str): Color map to use for the heatmap (default is "viridis").
    
    Returns:
    - None: Displays the heatmap.
    """
    # Load the data
    genome_profile = pd.read_csv(file_path, sep='\t')

    # Set the gene names as columns and spots as rows
    genome_profile = genome_profile.set_index('Unnamed: 0').transpose()

    # Perform hierarchical clustering on the spots
    linkage_matrix = linkage(genome_profile, method=clustering_method)

    # Plot the heatmap with hierarchical clustering on the spots (rows)
    g = sns.clustermap(
        genome_profile,
        row_cluster=True,
        col_cluster=False,
        row_linkage=linkage_matrix,
        cmap=cmap,
        figsize=(15, 10)
    )

    # Adjust the layout to make room for the title
    plt.subplots_adjust(top=0.92)

    # Set the centered title
    g.fig.suptitle(
        "Genome Profile Heatmap with Hierarchical Clustering on Spots",
        x=0.5,
        y=0.98,
        ha='center',
        fontsize=16
    )

    plt.show()

plot_genome_profile_heatmap("/Users/patricktruong/git/st_phylo/data/synthetic/spatial_data/st-genome_profile.tsv")

def plot_spatial_with_histology(file_paths, histology_image_path, gene_name, spot_size=3, scalef=3.4, dx=52, dy=45, alpha_img=0.7):
    """
    Plots annotation, gene expression, and cell count on a histology image background.

    Parameters:
    - file_paths (dict): Paths to the metadata, expression, and annotation TSV files as a dictionary.
                         Example: {'meta': 'path/to/meta.tsv', 'exp': 'path/to/expression.tsv', 'annotation': 'path/to/annotation.tsv'}
    - histology_image_path (str): Path to the histology image file.
    - gene_name (str): Gene name to display in the center plot for gene expression.
    - spot_size (int): Size of spots in the plot (default is 3).
    - scalef (float): Scale factor for the image (default is 3.4).
    - dx (int): Horizontal shift for alignment (default is 52).
    - dy (int): Vertical shift for alignment (default is 45).
    - alpha_img (float): Transparency of the background image (default is 0.7).
    """

    # Load metadata, expression data, and annotations
    df = pd.read_csv(file_paths['meta'], sep="\t")
    exp = pd.read_csv(file_paths['exp'], sep="\t")
    annotations = pd.read_csv(file_paths['annotation'], sep="\t", header=None, names=["Spot", "annotation"])

    # Process expression data
    exp.set_index('Unnamed: 0', inplace=True)
    expression_data = exp.T

    # Create AnnData object with expression data
    adata = AnnData(X=expression_data.values)
    adata.var_names = expression_data.columns  # Set gene names in var_names
    adata.obs['cell_count'] = df['cell_count'].values
    adata.obsm['spatial'] = df[['x', 'y']].values
    adata.obs['annotation'] = annotations['annotation'].values  # Add annotation data

    # Load histology image and add it to adata
    img = Image.open(histology_image_path)
    img = img.transpose(Image.FLIP_TOP_BOTTOM)

    # Set up spatial information in adata
    adata.uns["spatial"] = {
        "synthetic_dataset": {
            "images": {"hires": np.array(img)},
            "scalefactors": {
                "spot_diameter_fullres": 0.1,
                "tissue_hires_scalef": scalef,
            }
        }
    }

    # Adjust spot coordinates
    adata.obsm['spatial'][:, 0] += dx
    adata.obsm['spatial'][:, 1] += dy

    # Create a grid plot
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))

    # Plot annotation on the left
    sc.pl.spatial(
        adata,
        color='annotation',
        spot_size=spot_size,
        img_key="hires",
        library_id="synthetic_dataset",
        alpha_img=alpha_img,
        ax=axes[0],
        show=False
    )
    axes[0].set_title("Annotation")

    # Plot gene expression in the middle
    sc.pl.spatial(
        adata,
        color=gene_name,
        spot_size=spot_size,
        img_key="hires",
        library_id="synthetic_dataset",
        alpha_img=alpha_img,
        ax=axes[1],
        show=False
    )
    axes[1].set_title(f"Gene Expression: {gene_name}")

    # Plot cell count on the right
    sc.pl.spatial(
        adata,
        color='cell_count',
        spot_size=spot_size,
        img_key="hires",
        library_id="synthetic_dataset",
        alpha_img=alpha_img,
        ax=axes[2],
        show=False
    )
    axes[2].set_title("Cell Count")

    plt.tight_layout()
    plt.show()


plot_spatial_with_histology(file_paths, histology_image_path, gene_name, spot_size=3, scalef=3.4, dx=52, dy=45, alpha_img=0.7)
    # Load metadata, expression data, and annotations
df = pd.read_csv("/Users/patricktruong/git/st_phylo/data/synthetic/spatial_data/st-meta.tsv", sep="\t")
exp = pd.read_csv("/Users/patricktruong/git/st_phylo/data/synthetic/spatial_data/st-expression.tsv", sep="\t")
annotations = pd.read_csv("/Users/patricktruong/git/st_phylo/data/synthetic/spatial_data/st-annotation.tsv", sep="\t", header=None, names=["Spot", "annotation"])


# Define base paths for data and images
base_path = "/Users/patricktruong/git/st_phylo/data/synthetic/spatial_data"
img_path = "/Users/patricktruong/git/st_phylo/data/synthetic/visual"

# Adjust file paths using base_path and img_path
file_paths = {
    'meta': f"{base_path}/st-meta.tsv",
    'exp': f"{base_path}/st-expression.tsv",
    'annotation': f"{base_path}/st-annotation.tsv"
}
histology_image_path = f"{img_path}/tissue.png"

# Call the function with the updated paths
plot_spatial_with_histology(file_paths, histology_image_path, gene_name="Gene_0")
