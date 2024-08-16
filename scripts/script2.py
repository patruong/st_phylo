import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc
import squidpy as sq


adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

adata
adata2

adata.uns["spatial"]["V1_Human_Lymph_Node"].keys()
adata.uns["spatial"]["V1_Human_Lymph_Node"]["scalefactors"]
adata.uns["spatial"]["V1_Human_Lymph_Node"]["metadata"]

# read in from scipt.py
adata2.uns["spatial"]["Patient_2_H2_2"].keys()
adata2.uns["spatial"]["Patient_2_H2_2"]["scalefactors"]
