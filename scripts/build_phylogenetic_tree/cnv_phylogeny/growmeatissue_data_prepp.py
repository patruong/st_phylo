import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec  # Import gridspec

from scipy.cluster.hierarchy import linkage, dendrogram
import scanpy as sc
import numpy as np
from anndata import AnnData
import squidpy as sq
from PIL import Image



def plot_genome_profile_heatmaps_with_dendrograms(file_path, clustering_method='ward', cmap="viridis"):
    """
    Plots two genome profile heatmaps side by side, each with their dendrograms and spot annotations:
    - Left: Without log transformation.
    - Right: With log-transformed gene expression.

    Parameters:
    - file_path (str): Path to the genome profile TSV file.
    - clustering_method (str): Clustering method to use for hierarchical clustering (default is 'ward').
    - cmap (str): Color map to use for the heatmaps (default is "viridis").

    Returns:
    - None: Displays the heatmaps with dendrograms.
    """
    # Load the data
    genome_profile = pd.read_csv(file_path, sep='\t')

    # Set the gene names as columns and spots as rows
    genome_profile = genome_profile.set_index(genome_profile.columns[0]).transpose()

    # Filter away spots where all genes are 1.0
    genome_profile = genome_profile.loc[~(genome_profile.iloc[:, 1:] == 1.0).all(axis=1)]

    # Apply log transformation (adding 1 to avoid log(0))
    genome_profile_log = np.log1p(genome_profile)

    # Prepare the datasets in a list
    datasets = [
        ('Without Log Transformation', genome_profile),
        ('With Log Transformation', genome_profile_log)
    ]

    # Create a figure
    fig = plt.figure(figsize=(20, 10))
    fig.suptitle(
        "Genome Profile Heatmaps with Hierarchical Clustering on Spots",
        fontsize=20,
        y=0.98  # Adjust the y position if necessary
    )

    # Create gridspec
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])

    for i, (title, data) in enumerate(datasets):
        # Perform hierarchical clustering on the spots (rows)
        linkage_matrix = linkage(data, method=clustering_method)

        # Create a gridspec within the main gridspec for dendrogram and heatmap
        gs_inner = gridspec.GridSpecFromSubplotSpec(
            1, 2,
            width_ratios=[0.3, 1],
            wspace=0.05,  # Adjust spacing between dendrogram and heatmap
            subplot_spec=gs[i]
        )

        # Plot dendrogram
        ax_dendro = fig.add_subplot(gs_inner[0])
        dendro = dendrogram(
            linkage_matrix,
            orientation='left',
            no_labels=True,
            ax=ax_dendro,
            color_threshold=0
        )
        ax_dendro.invert_yaxis()
        ax_dendro.axis('off')

        # Get the order of the leaves
        row_order = dendro['leaves']

        # Reorder the data according to clustering
        data_ordered = data.iloc[row_order, :]

        # Plot heatmap
        ax_heatmap = fig.add_subplot(gs_inner[1])
        sns.heatmap(
            data_ordered,
            ax=ax_heatmap,
            cmap=cmap,
            cbar_kws={
                'label': (
                    'Expression Level' if title == 'Without Log Transformation'
                    else 'Log(Expression Level)'
                )
            },
            xticklabels=True,
            yticklabels=True  # Display spot annotations
        )
        ax_heatmap.set_title(title, fontsize=16)
        ax_heatmap.set_xlabel("Genes")
        ax_heatmap.set_ylabel("")

        # Adjust the heatmap's y-axis to match the dendrogram
        ax_heatmap.yaxis.set_ticks_position('right')
        ax_heatmap.yaxis.set_label_position('right')
        ax_heatmap.tick_params(axis='y', which='both', length=0)
        # Set the y-tick labels to the spot names
        ax_heatmap.set_yticklabels(data_ordered.index.tolist(), fontsize=8)
        # Optionally adjust the rotation
        plt.setp(ax_heatmap.get_yticklabels(), rotation=0)

    # Adjust layout to make room for the title
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

# FILTER AWAY CNV WITHOUT ANY CHANGE...

file_path = "/Users/patricktruong/git/st_phylo/data/synthetic/spatial_data/st-genome_profile.tsv"

genome_profile = pd.read_csv(file_path, sep='\t')
genome_profile = genome_profile.set_index(genome_profile.columns[0]).transpose()

plot_genome_profile_heatmaps_with_dendrograms("/Users/patricktruong/git/st_phylo/data/synthetic/spatial_data/st-genome_profile.tsv")



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

def plot_spatial_with_histology(file_paths, histology_image_path, gene_name, spot_size=3, scalef=3.4, dx=52, dy=45, alpha_img=0.7, reannotate_aberrant = True):
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

    # Reannotate
    if reannotate_aberrant == True:
        annotations = reannotate(annotations, df)
        annotations["annotation"] = annotations['new_annotation']

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
plot_spatial_with_histology(file_paths, histology_image_path, gene_name="Gene_0", reannotate_aberrant = True)
plot_spatial_with_histology(file_paths, histology_image_path, gene_name="Gene_0", reannotate_aberrant = False)


### Map spot to aberrant group 

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

annotations = reannotate(annotations, df)
annotations.new_annotation.unique()

annotations["annotation"] = annotations['new_annotation']


annotations = reannotate(annotations, df)
annotations["annotation"] = annotations['new_annotation']
