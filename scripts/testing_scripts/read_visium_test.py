
import scanpy as sc

path = "/home/ptruong/git/st_phylo/test"

adata = sp.read_visium(path)
adata.var_names_make_unique()
adata