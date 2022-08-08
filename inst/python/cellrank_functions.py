# functions to fit transition matrices using cellrank
from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel, PseudotimeKernel
from cellrank.tl.estimators import GPCCA, CFLARE

def fit_cellrank_velocity_kernel(adata, velocity_weight = 0.8):
    """
    fit fine-grained (i.e. cell-to-cell) transition matrix based on 
    RNA velocity estimates
    
    :args adata: AnnData object with gene counts and unspliced/spliced counts added
    :args velocity_weight: float, weight given to velocity-kernel transition matrix. between 0 and 1. (1 - velocity_weight) will be the weight given to the kNN graph constructed based on the GEX data.
    :returns: a dictionary with two items:
      - key 'cellrank_obj': cellrank.tl.estimators.CFLARE object. with attribute 'transition_matrix' which gives the cell-to-cell transition matrix of class sp.sparse.csr_matrix.
      - key 'transition_matrix': dense transition matrix inside `cellrank_obj` for convenience.      - key 'CellID': list, cell identifiers for each row/column of `transition_matrix`.
    """
    vk = VelocityKernel(adata).compute_transition_matrix()
    ck = ConnectivityKernel(adata).compute_transition_matrix()
    # combine kernel:
    combined_kernel = velocity_weight * vk + (1 - velocity_weight) * ck
    # fit markov chain defined by this combined_kernel transition matrix
    g = CFLARE(combined_kernel)
    g.fit()
    return {'cellrank_obj': g, 'transition_matrix': g.transition_matrix.todense(), 'CellID': [str(i) for i in adata.obs_names]}

def fit_cellrank_pseudotime_kernel(adata, pseudotime_key, pseudotime_weight = 0.8):
    """
    fit fine-grained (i.e. cell-to-cell) transition matrix based on 
    a pseudotime ordering
    
    :args adata: AnnData object with GEX data. pseudotime_key is found in adata.obs.
    :args pseudotime_key: str, name of column in adata.obs which holds the pseudotime ordering.
    :args pseudotime_weight: float, weight given to pseudotime-kernel transition matrix. between 0 and 1. (1 - pseudotime_weight) will be the weight given to the kNN graph constructed based on the GEX data.
    
    :returns: a dictionary with two items:
      - key 'cellrank_obj': cellrank.tl.estimators.CFLARE object. with attribute 'transition_matrix' which gives the cell-to-cell transition matrix of class sp.sparse.csr_matrix.
      - key 'transition_matrix': dense transition matrix inside `cellrank_obj` for convenience.
      - key 'CellID': list, cell identifiers for each row/column of `transition_matrix`.
    """
    vk = PseudotimeKernel(adata, time_key = pseudotime_key).compute_transition_matrix()
    ck = ConnectivityKernel(adata).compute_transition_matrix()
    # combine kernel:
    combined_kernel = pseudotime_weight * vk + (1 - pseudotime_weight) * ck
    # fit markov chain defined by this combined_kernel transition matrix
    g = CFLARE(combined_kernel)
    g.fit()
    return {'cellrank_obj': g, 'transition_matrix': g.transition_matrix.todense(), 'CellID': [str(i) for i in adata.obs_names]}



