import pandas as pd

df = pd.read_csv("Synthethic_data/spatial_data/st-expression.tsv", sep = "\t", index_col = 0)
df = pd.read_csv("Synthethic_data/spatial_data/st-annotation.tsv", sep = "\t", index_col = 0)
df = pd.read_csv("Synthethic_data/spatial_data/st-genome_profile.tsv", sep = "\t", index_col = 0)

# need to have inferCNV data on this.
