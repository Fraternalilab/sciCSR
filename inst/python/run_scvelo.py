import scanpy as sc
import scvelo as scv
import argparse
import sys
scv.settings.verbosity = 3

if __name__=="__main__":
        parser = argparse.ArgumentParser(description="run scVelo")
        parser.add_argument('--anndata_file', help = 'input filepath to AnnData .h5ad file', type = str, required = True)
        parser.add_argument('--anndata_out_file', help = 'output filepath to AnnData .h5ad file', type = str, required = True)
        parser.add_argument('--mode', help = 'mode', default = 'dynamical', choices = ['deterministic', 'stochastic', 'dynamical'])
        parser.add_argument('--min_shared_counts', help = 'Minimum number of counts (both unspliced and spliced) required for a gene.', default = int(20))
        parser.add_argument('--n_top_genes', help = 'Number of genes to keep', default = int(2000))
        parser.add_argument('--reduction', help = 'dimensionality reduction to use', type = str, default = 'umap')
        args = parser.parse_args()
        
        min_shared_counts = int(args.min_shared_counts)
        n_top_genes = int(args.n_top_genes)

        adata = sc.read_h5ad(args.anndata_file)
        print("running scVelo pipeline ...")
        # first calculate moments
        scv.pp.filter_and_normalize(adata, min_shared_counts = min_shared_counts, n_top_genes = n_top_genes)
        scv.pp.neighbors(adata)
        scv.pp.moments(adata, n_pcs = None, n_neighbors = None)
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode = args.mode)
        scv.tl.velocity_graph(adata)
        scv.tl.velocity_embedding(adata, basis = args.reduction, autoscale = False)
        adata.write_h5ad(args.anndata_out_file)
        
