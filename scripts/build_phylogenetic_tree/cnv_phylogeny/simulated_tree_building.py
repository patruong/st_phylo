import pandas as pd
from anytree import Node, RenderTree
from anytree.exporter import DotExporter
import matplotlib.pyplot as plt
from PIL import Image
import subprocess
import tempfile


# Load and remove non-CNV altered spots (not interesting)

file_path = "/Users/patricktruong/git/st_phylo/data/synthetic/spatial_data/st-genome_profile.tsv"

genome_profile = pd.read_csv(file_path, sep='\t')
genome_profile = genome_profile.set_index(genome_profile.columns[0]).transpose()
genome_profile = genome_profile.loc[~(genome_profile.iloc[:, 1:] == 1.0).all(axis=1)]

genome_profile.index

# Step 1: Load and Parse the Lineage Data
lineage_df = pd.read_csv('/Users/patricktruong/git/st_phylo/data/synthetic/lineage.tsv', sep='\t')
lineage_df = lineage_df.rename(columns={"child": "Child", "parent": "Parent"}).loc[:, ["Child", "Parent"]]
lineage_df['Child'] = lineage_df['Child'].astype(int)
lineage_df['Parent'] = lineage_df['Parent'].astype(int)
parent_child_map = lineage_df.set_index('Child')['Parent'].to_dict()

# Step 2: Load and Parse the Spot Data
spots_df = pd.read_csv('/Users/patricktruong/git/st_phylo/data/synthetic/spatial_data/cell_by_spot.txt', sep=' : ', 
                       engine='python', header=None, names=["Spot", "Cells"])
spots_df['Spot'] = spots_df['Spot'].str.strip()
spots_df['Cells'] = spots_df['Cells'].apply(lambda x: list(map(int, x.split(','))))
spot_to_cells = dict(zip(spots_df['Spot'], spots_df['Cells']))

spots_df = spots_df[spots_df['Spot'].isin(genome_profile.index)]

spots_df


