import scanpy as sc
import scvelo as scv
import argparse
import sys
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel, PseudotimeKernel
scv.settings.figdir = ""
scv.settings.plot_prefix = ""

if __name__=="__main__":
        parser = argparse.ArgumentParser(description="plot a velocity plot using scvelo and save as png/svg/pdf")
        parser.add_argument('--anndata', help = 'filepath to AnnData .h5ad file', type = str, required = True)
        parser.add_argument('--type', help = 'use velocity/csr/shm information to plot arrows.', type = str, choices = ['velocity', 'shm', 'csr'], default = 'velocity')
        parser.add_argument('--color', help = 'Key for annotations of observations/cells or variables/genes', default = None)
        parser.add_argument('--palette', help = 'comma-separated list of colours to use in the plot. Length must correspond to the number of levels in the argument given by \'color\'.', type = str, default = None)
        parser.add_argument('--title', help = 'title of the plot', type = str, default = None)
        parser.add_argument('--output', help = 'output png/svg/pdf filename', type = str, required = True, default = 'velocity.png')
        parser.add_argument('--style', help = 'output arrows on grid (default) or as streams.', type = str, default = 'grid', choices = ['grid', 'stream'])
        parser.add_argument('--dpi', help = 'dots per inch/ Ignored if outputting svg/pdf.', type = int, default = 600)
        parser.add_argument('--components', help = 'component of dim reduction to show. e.g. put \'1,3\' if desired UMAP plot displays the UMAP_1 and UMAP_3 axes.', type = str, default = '1,2')
        args = parser.parse_args()

        adata = sc.read_h5ad(args.anndata)

        velocity_key = 'velocity'
        if args.type != 'velocity':
            if args.type == 'csr':
                vk = PseudotimeKernel(adata, time_key = 'csr_pot').compute_transition_matrix()
            elif args.type == 'shm':
                vk = PseudotimeKernel(adata, time_key = 'shm').compute_transition_matrix()
            vk.compute_projection()    
            velocity_key = 'T_fwd'

        if args.palette is not None:
            colour_list = args.palette.split(',')
        else:
            colour_list = None
        
        if args.title is not None:
            main = args.title
        else:
            main = args.type

        if args.style == 'grid':
                scv.pl.velocity_embedding_grid(adata, vkey = velocity_key, title = main, save = args.output, dpi = args.dpi, color = args.color, palette = colour_list, components = args.components, arrow_length=3, scale = 2, show = False)
        elif args.style == 'stream':
                scv.pl.velocity_embedding_stream(adata, vkey = velocity_key, title = main, save = args.output, dpi = args.dpi, color = args.color, palette = colour_list, components = args.components, show = False)

