#!/usr/bin/env python3

import argparse
import sys

# %%
import scanpy as sc
import plotnine as p9

import liana as li
import cell2cell as c2c
import decoupler as dc # needed for pathway analysis

import warnings
warnings.filterwarnings('ignore')
from collections import defaultdict

# %%
import pandas as pd
import matplotlib.pyplot as plt
# %matplotlib inline
import plotnine as p9
import seaborn as sns

def parse_args(args=None):
    Description = "Liana cell2cell"
    Epilog = "Example usage: python 001.liana.py <adata> <adata_out>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--adata")
    parser.add_argument("--adata_out")

    return parser.parse_args(args)

def liana(adata,adata_out):
# %%
    #adata = adata[adata.obs.sample_type.isin(["tumor","normal"])] #subset tumor and normal
    adata.layers["log1p_norm"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e6, layer="log1p_norm")
    sc.pp.log1p(adata, base=6, layer="log1p_norm")

    # %% [markdown]
    # ## LIANA

    # %%
    # select ligand receptor resources
    # LIANA uses human gene symbol identifiers,

    # %%
    #lr_pairs = li.resource.select_resource("consensus")

    # %%
   #li.mt.rank_aggregate.by_sample(adata, 
   #                            sample_key = "sample_id",
   #                            groupby = "cell_type",
   #                            resource_name ="consensus",
   #                            expr_prop = 0.1,
   #                            min_cells = 5,
   #                            n_perms = 100,
   #                            use_raw=False,
   #                            layer = "log1p_norm",
   #                            verbose = True, 
   #                            inplace = True)
    
    li.mt.rank_aggregate(adata,
                          groupby='cell_type_fine',
                            expr_prop=0.1,
                            resource_name='consensus',
                          verbose=True,
                          key_added='rank_aggregate',
                            layer = "log1p_norm",
                            n_jobs=7,
                              use_raw = False)

    # %%
    adata.uns['rank_aggregate'].to_csv("LIANA_by_sample.csv", index = False)
    adata.write_h5ad(adata_out, compression="gzip")
    adata_out=adata

    return adata_out



def main(args=None):
    args = parse_args(args)

    liana(
        adata=sc.read_h5ad(args.adata),
        adata_out=args.adata_out,
    )
   

if __name__ == "__main__":
    sys.exit(main())