#!/usr/bin/env python3

import argparse
import sys
import scanpy as sc


import liana as li
import cell2cell as c2c


import warnings
warnings.filterwarnings('ignore')
from collections import defaultdict




def parse_args(args=None):
    Description = "Liana cell2cell"
    Epilog = "Example usage: python 002.tensor.py <adata> <adata_out> <tensor_out> <tensor_meta_out> "

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--adata")
    parser.add_argument("--build_tensor_out")
    parser.add_argument("--meta_build_tensor_out")
    parser.add_argument("--run_tensor_out")
    parser.add_argument("--outdir")
    return parser.parse_args(args)


def tensor_cell2cell(adata, build_tensor_out, meta_build_tensor_out,run_tensor_out, outdir):
    ######################### tensor_cell2cell
   
    adata = sc.read_h5ad(adata)

    ## LIANA BUILD TENSOR 
    sample_key = 'sample_id'
    condition_key = 'sample_type'
    groupby = 'cell_type'

    build_tensor = li.multi.to_tensor_c2c(adata,
                                    sample_key=sample_key,
                                    score_key='magnitude_rank', # can be any score from liana
                                    how='outer_cells' # how to join the samples
                                    )
    
    print(build_tensor.shape)
    #export tensor and meta tensor 
    c2c.io.save_data.export_variable_with_pickle(variable=build_tensor, filename=build_tensor_out)

    context_dict = adata.obs[[sample_key, condition_key]].drop_duplicates()
    context_dict = dict(zip(context_dict[sample_key], context_dict[condition_key]))
    context_dict = defaultdict(lambda: 'Unknown', context_dict)
  

    meta_build_tensor = c2c.tensor.generate_tensor_metadata(interaction_tensor=build_tensor,
                                                  metadata_dicts=[context_dict, None, None, None],
                                                  fill_with_order_elements=True
                                                  )
    
     
    c2c.io.save_data.export_variable_with_pickle(variable=meta_build_tensor, filename=meta_build_tensor_out)

    run_tensor = c2c.analysis.run_tensor_cell2cell_pipeline(build_tensor,
                                                    meta_build_tensor,
                                                    copy_tensor=True, # Whether to output a new tensor or modifying the original
                                                    rank=6, # Number of factors to perform the factorization. If None, it is automatically determined by an elbow analysis. Here, it was precomuputed.
                                                    tf_optimization='robust', # To define how robust we want the analysis to be.
                                                    random_state=0, # Random seed for reproducibility
                                                    device='cpu', # Device to use. If using GPU and PyTorch, use 'cuda'. For CPU use 'cpu'
                                                    elbow_metric='error', # Metric to use in the elbow analysis.
                                                    smooth_elbow=False, # Whether smoothing the metric of the elbow analysis.
                                                    upper_rank=20, # Max number of factors to try in the elbow analysis
                                                    tf_init='random', # Initialization method of the tensor factorization
                                                    tf_svd='numpy_svd', # Type of SVD to use if the initialization is 'svd'
                                                    cmaps=None, # Color palettes to use in color each of the dimensions. Must be a list of palettes.
                                                    sample_col='Element', # Columns containing the elements in the tensor metadata
                                                    group_col='Category', # Columns containing the major groups in the tensor metadata
                                                    output_fig=False, # Whether to output the figures. If False, figures won't be saved a files if a folder was passed in output_folder.
                                                    )
 

    c2c.io.save_data.export_variable_with_picke(variable=run_tensor, filename=run_tensor_out)

    return build_tensor, meta_build_tensor, run_tensor
    


def main(args=None):
    args = parse_args(args)

    tensor_cell2cell(adata = args.adata, build_tensor_out = args.build_tensor_out, meta_build_tensor_out = args.meta_build_tensor_out, run_tensor_out = args.run_tensor_out,outdir = args.outdir)


if __name__ == "__main__":
    sys.exit(main())