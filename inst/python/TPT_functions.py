# customs functions to override committor (fwd/bwd) calculation in deeptime
import multiprocess
from functools import partial
import deeptime.markov.tools.analysis as msmana
from deeptime.markov.msm import MarkovStateModel
from deeptime.markov.tools.analysis._api import is_reversible
import deeptime.markov.tools.flux as tptapi
from deeptime.markov._reactive_flux import ReactiveFlux
from deeptime.util.types import ensure_array, ensure_number_array, ensure_integer_array, ensure_floating_array
import numpy as np
from scipy.linalg import solve
from scipy.sparse.linalg import spsolve
import scipy.sparse as sparse
from scipy.stats import norm
import os
import random
import pandas as pd

def committor(T, A, B, forward=True, mu=None):
    r"""Compute the committor between sets of microstates.
    The committor assigns to each microstate a probability that being
    at this state, the set B will be hit next, rather than set A
    (forward committor), or that the set A has been hit previously
    rather than set B (backward committor). See :footcite:`metzner2009transition` for a
    detailed mathematical description. The present implementation
    uses the equations given in :footcite:`noe2009constructing`.
    Parameters
    ----------
    T : (M, M) ndarray or scipy.sparse matrix
        Transition matrix
    A : array_like
        List of integer state labels for set A
    B : array_like
        List of integer state labels for set B
    forward : bool
        If True compute the forward committor, else
        compute the backward committor.
    Returns
    -------
    q : (M,) ndarray
        Vector of comittor probabilities.
    Notes
    -----
    Committor functions are used to characterize microstates in terms
    of their probability to being visited during a reaction/transition
    between two disjoint regions of state space A, B.
    **Forward committor**
    The forward committor :math:`q^{(+)}_i` is defined as the probability
    that the process starting in `i` will reach `B` first, rather than `A`.
    Using the first hitting time of a set :math:`S`,
    .. math:: T_{S}=\inf\{t \geq 0 | X_t \in S \}
    the forward committor :math:`q^{(+)}_i` can be fromally defined as
    .. math:: q^{(+)}_i=\mathbb{P}_{i}(T_{A}<T_{B}).
    The forward committor solves to the following boundary value problem
    .. math::  \begin{array}{rl} \sum_j L_{ij} q^{(+)}_{j}=0 & i \in X \setminus{(A \cup B)} \\
                q_{i}^{(+)}=0 & i \in A \\
                q_{i}^{(+)}=1 & i \in B
                \end{array}
    :math:`L=T-I` denotes the generator matrix.
    **Backward committor**
    The backward committor is defined as the probability that the process
    starting in :math:`x` came from :math:`A` rather than from :math:`B`.
    Using the last exit time of a set :math:`S`,
    .. math:: t_{S}=\sup\{t \geq 0 | X_t \notin S \}
    the backward committor can be formally defined as
    .. math:: q^{(-)}_i=\mathbb{P}_{i}(t_{A}<t_{B}).
    The backward comittor solves another boundary value problem
    .. math::  \begin{array}{rl}
                \sum_j K_{ij} q^{(-)}_{j}=0 & i \in X \setminus{(A \cup B)} \\
                q_{i}^{(-)}=1 & i \in A \\
                q_{i}^{(-)}=0 & i \in B
                \end{array}
    :math:`K=(D_{\pi}L)^{T}` denotes the adjoint generator matrix.
    References
    ----------
    .. footbibliography::
    Examples
    --------
    >>> import numpy as np
    >>> from deeptime.markov.tools.analysis import committor
    >>> T = np.array([[0.89, 0.1, 0.01], [0.5, 0.0, 0.5], [0.0, 0.1, 0.9]])
    >>> A = [0]
    >>> B = [2]
    >>> u_plus = committor(T, A, B)
    >>> u_plus
    array([0. , 0.5, 1. ])
    >>> u_minus = committor(T, A, B, forward=False)
    >>> u_minus
    array([1.        , 0.45454545, 0.        ])
    """
    T = ensure_number_array(T, ndim=2)
    A = ensure_integer_array(A, ndim=1)
    B = ensure_integer_array(B, ndim=1)
    if forward:
        return forward_committor(T, A, B)
    else:
        """ if P is time reversible backward commitor is equal 1 - q+"""
        if is_reversible(T, mu=mu):
            return 1.0 - forward_committor(T, A, B)

        else:
            return backward_committor(T, A, B)
        
def _set_up_linear_system(K, A, B):
    """Assemble left-hand side W for linear system"""
    import scipy.sparse as sparse
    import numpy as np
    """Equation (I)"""
    W = 1.0 * K
    """Equation (II)"""
    if sparse.issparse(W):
        W[list(A), :] = 0.0
        W = W + sparse.coo_matrix((np.ones(len(A)), (list(A), list(A))), shape=W.shape).tocsr()
    else:
        W[list(A), :] = 0.0
        W[list(A), list(A)] = 1.0
    """Equation (III)"""
    if sparse.issparse(W):
        W[list(B), :] = 0.0
        W = W + sparse.coo_matrix((np.ones(len(B)), (list(B), list(B))), shape=W.shape).tocsr()
    else:
        W[list(B), :] = 0.0
        W[list(B), list(B)] = 1.0
    return W


def forward_committor(T, A, B):
    r"""Forward committor between given sets.
    The forward committor u(x) between sets A and B is the probability
    for the chain starting in x to reach B before reaching A.
    Parameters
    ----------
    T : (M, M) ndarray
        Transition matrix
    A : array_like
        List of integer state labels for set A
    B : array_like
        List of integer state labels for set B
    Returns
    -------
    u : (M, ) ndarray
        Vector of forward committor probabilities
    Notes
    -----
    The forward committor is a solution to the following
    boundary-value problem
    .. math::
        \sum_j L_{ij} u_{j}=0    for i in X\(A u B) (I)
                      u_{i}=0    for i \in A        (II)
                      u_{i}=1    for i \in B        (III)
    with generator matrix L=(P-I).
    """
    import numpy as np
    from scipy.linalg import solve
    from scipy.sparse.linalg import spsolve
    import scipy.sparse as sparse

    A = set(A)
    B = set(B)
    AB = A.intersection(B)
    if len(AB) > 0:
        raise ValueError("Sets A and B have to be disjoint")
    L = T - np.eye(T.shape[0])  # Generator matrix

    W = _set_up_linear_system(L, A, B)
    """Assemble right hand side r for linear system"""
    """Equation (I+II)"""
    r = np.zeros(T.shape[0])
    """Equation (III)"""
    r[list(B)] = 1.0

    u = solve(W, r) if not sparse.issparse(W) else spsolve(W, r)
    return u


def backward_committor(T, A, B, mu=None):
    r"""Backward committor between given sets.
    The backward committor u(x) between sets A and B is the
    probability for the chain starting in x to have come from A last
    rather than from B.
    Parameters
    ----------
    T : (M, M) ndarray
        Transition matrix
    A : array_like
        List of integer state labels for set A
    B : array_like
        List of integer state labels for set B
    mu : (M, ) ndarray (optional)
        Stationary vector
    Returns
    -------
    u : (M, ) ndarray
        Vector of forward committor probabilities
    Notes
    -----
    The forward committor is a solution to the following
    boundary-value problem
    .. math::
        \sum_j K_{ij} \pi_{j} u_{j}=0    for i in X\(A u B) (I)
                                  u_{i}=1    for i \in A        (II)
                                  u_{i}=0    for i \in B        (III)
    with adjoint of the generator matrix K=(D_pi(P-I))'.
    """
    import deeptime.markov.tools.analysis as msmana
    import numpy as np
    from scipy.linalg import solve
    from scipy.sparse.linalg import spsolve
    import scipy.sparse as sparse

    A = set(A)
    B = set(B)
    AB = A.intersection(B)
    if len(AB) > 0:
        raise ValueError("Sets A and B have to be disjoint")
    if mu is None:
        mu = msmana.stationary_distribution(T)
    if sparse.issparse(T):
        L = T - sparse.eye(T.shape[0], T.shape[0])
        D = sparse.diags([mu, ], [0, ])
        K = (D.dot(L)).T
    else:
        K = np.transpose(mu[:, np.newaxis] * (T - np.eye(T.shape[0])))

    """Assemble left-hand side W for linear system"""
    W = _set_up_linear_system(K, A, B)
    """Assemble right-hand side r for linear system"""
    """Equation (I)+(III)"""
    r = np.zeros(T.shape[0])
    """Equation (II)"""
    r[list(A)] = 1.0

    u = solve(W, r) if not sparse.issparse(W) else spsolve(W, r)

    return u        

def fit_tpt(transition_matrix, source_states, target_states, verbose = True):
    """
    Rewrote from deeptime.markov._reactive_flux.reactive_flux
    slimming down memory usage of some linear algebra steps

    """
    import deeptime.markov.tools.analysis as msmana
    from deeptime.markov.msm import MarkovStateModel
    from deeptime.markov.tools.analysis._api import is_reversible
    import deeptime.markov.tools.flux as tptapi
    from deeptime.markov._reactive_flux import ReactiveFlux
    from deeptime.util.types import ensure_array, ensure_number_array, ensure_integer_array, ensure_floating_array
    import numpy as np
    from scipy.linalg import solve
    from scipy.sparse.linalg import spsolve
    import scipy.sparse as sparse
    from scipy.stats import norm
    source_states = ensure_array(source_states, dtype=int)
    target_states = ensure_array(target_states, dtype=int)

    if len(source_states) == 0 or len(target_states) == 0:
        raise ValueError('set A or B is empty')

    n_states = transition_matrix.shape[0]
    if len(source_states) > n_states or len(target_states) > n_states \
            or max(source_states) > n_states or max(target_states) > n_states:
        raise ValueError('set A or B defines more states than the given transition matrix.')
        
    if verbose:
        print("Input check complete. Now compute ...")

    # stationary distribution    
    if verbose:
        print("Stationary distribution ...")
    stationary_distribution = msmana.stationary_distribution(transition_matrix)
    
    # forward committor
    if verbose:
        print("forward committor ...")
    qplus = forward_committor(transition_matrix, source_states, target_states)
    
    # backward committor
    if verbose:
        print("backward committor ...")
    if msmana.is_reversible(transition_matrix, mu=stationary_distribution):
        qminus = 1.0 - qplus
    else:
        qminus= backward_committor(transition_matrix, source_states, target_states,
                                   mu = stationary_distribution)
    
    # gross flux
    if verbose:
        print("gross flux ...")
    grossflux = tptapi.flux_matrix(transition_matrix, stationary_distribution, 
                                   qminus, qplus, netflux=False)
    
    # net flux
    if verbose:
        print("net flux ...")
    netflux = tptapi.to_netflux(grossflux)
    
    return ReactiveFlux(source_states, target_states, net_flux=netflux,
                        stationary_distribution=stationary_distribution,
                        qminus=qminus, qplus=qplus, gross_flux=grossflux)

def getReshuffledFlux_(i, tm, cluster_ident, source, target):
    import random
    import numpy as np
    random.seed(i + 1)
    # random reshuffle columns
    indices = np.arange(tm.shape[1])
    random.shuffle(indices)
    # tm = tm[:, indices]
    cluster_ident_random = dict()
    j = 0
    for key, item in cluster_ident.items():
        cluster_ident_random[key] = indices[j:(j + len(item))]
        j += len(item)
    flux = fit_tpt(tm, cluster_ident_random[source], cluster_ident_random[target], verbose=False)
    (tpt_sets, tpt_coarse) = flux.coarse_grain([np.array(item) for c, item in cluster_ident_random.items()])
    # flux = fit_tpt(tm, cluster_ident[source], cluster_ident[target], verbose=False)
    # (tpt_sets, tpt_coarse) = flux.coarse_grain([np.array(item) for c, item in cluster_ident.items()])
    # get the gross_flux matrix
    gross_flux = tpt_coarse.gross_flux
    gross_flux = gross_flux * 100 / gross_flux.sum()
    # get the net_flux matrix
    net_flux = tpt_coarse.net_flux
    net_flux = net_flux * 100 / net_flux.sum()
    return {"tpt": tpt_coarse, "gross_flux": gross_flux, "net_flux": net_flux}

def init_worker(sfit_tpt, s_set_up_linear_system, sforward_committor, sbackward_committor):
    """
    initializer individual worker in pool for getReshuffledFlux
    ref https://superfastpython.com/multiprocessing-pool-shared-global-variables/
    """
    import random
    import numpy as np
    global fit_tpt
    fit_tpt = sfit_tpt
    global _set_up_linear_system
    _set_up_linear_system = s_set_up_linear_system
    global forward_committor
    forward_committor = sforward_committor
    global backward_committor
    backward_committor = sbackward_committor

def getReshuffledFlux(transition_matrix, cluster_ident, source, target, n = 100, show_progress_bar = True):
    """
    reshuffle transition matrix columns to get 'randomised' flux
    
    :args trnasition_matrix: the transition matrix
    :args cluster_ident: the dictionary mapping the id of the transition matrix belonging to each cluster (ie macrostate)
    :args source: the key in cluster_ident representing the 'source' macrostate in TPT
    :args target: the key in cluster_ident representing the 'target' macrostate in TPT
    :args n: int, the number of random reshuffling
    
    :return: a list of length n, each a dictionary with 3 items:
    - key "tpt": coarse-grained TPT object
    - key "gross_flux": gross_flux from TPT, normalised as gross flux between i and j as % of flux across system
    - key "net_flux": net_flux from TPT, normalised as net flux between i and j as % of net flux across system
    """
    if show_progress_bar:
        try:
            import ipywidgets
            from tqdm.auto import tqdm
        except ImportError:
            try:
                from tqdm.std import tqdm
            except ImportError:
                tqdm = None
        
    out = []
    n_cpus = os.cpu_count()
    if n_cpus == 1:
        out = [getReshuffledFlux_(i, tm = transition_matrix, cluster_ident = cluster_ident, source = source, target = target) for i in range(int(n))]
    else:
        with multiprocess.Pool(int(n_cpus / 2), initializer = init_worker, initargs = (fit_tpt, _set_up_linear_system, forward_committor, backward_committor,)) as p:
            out = list(tqdm(p.imap(partial(getReshuffledFlux_, tm = transition_matrix, cluster_ident = cluster_ident, source = source, target = target), range(int(n))), total=int(n)))

    return out

def getStationaryDistribution(stationary_microstates, cluster_ident, seed):
    """
    calculate stationary distribution of macrostates (i.e. coarse-grained 
    Markov State Model) given a transition matrix.
    Used as function to compute stationary distribution of bootstrapped microstates for
    confidence intervals

    ;args stationary_microstates: (M, ) array, stationary distribution of the microstates.
    :args cluster_ident: dict, user-defined coarse-graining scheme. Keys are the cluster labels, and values are lists of row/col indices of transition_matrix which belongs to the given cluster label
    :args seed: int, random seed
    """
    random.seed(seed)
    out = dict()
    for key, item in cluster_ident.items():
        indices = item
        bs_indices = [random.choice(indices) for _ in indices]
        out[key] = np.sum([stationary_microstates[i] for i in bs_indices])
    return out    

def fit_coarse_grain_tpt(transition_matrix, cluster_ident, 
                         source_state, target_state, random_n = 100, verbose = True):
    """
    Applies Transition Path Theory (TPT) on a transition matrix and coarse grain 
    the Markov State Model (MSM) into macrostates under a user-defined scheme.
    
    :args transition_matrix: the cell-to-cell transition matrix
    :args cluster_ident: dict, user-defined coarse-graining scheme. Keys are the cluster labels, and values are lists of row/col indices of transition_matrix which belongs to the given cluster label
    :args source_state: str or list of str, key(s) in cluster_ident which correspond to the source state to fit TPT.
    :args target_state: str or list of str, key(s) in cluster_ident which correspond to the target state to fit TPT.
    :args random_n: int, number of times to reshuffle columns of transition_matrix to fit randomised TPT. For determining significance of observed flux between any two clusters. (Default: 100)
    :args verbose: Bool, should we print messages throughout? (Default: True)

    :returns: a dictionary with the following keys:
    - 'coarse_grain_tpt': coarse-grained ReactiveFlux object. TPT applied based on cluster_ident, source_state and target_state.
    - 'gross_flux': matrix of gross_flux from coarse_graind_tpt, normalised as % of total flux which flows from cluster i to j  
    - 'pathways': pd.DataFrame of possible pathways and their likelihood of traversal under the coarse-grained TPT, from source_state to target_state. A dataframe of two columns: 
      * 'path': pathway to flow from source_state to target_state, 
      * 'p_path': % of traversals following this pathway.
    - 'randomised_tpt': output of 'getReshuffledFlux'. coarse-grained TPT fitted based on random reshuffling of columns of transition_matrix.
    - 'significance': IOne-sided P-value that the observed gross_flux from cluster i to j is *greater* than those sampled in random_flux.
    - 'total_gross_flux': sum of all entries in `gross_flux`, but before the normalisation as percentages.
    - 'total_gross_flux_randomised': analogous to `total_gross_flux`, but for each iteration of `randomised_tpt`.
    - "bootstrap_stationary': stationary distribution in TPT models fitted on `random_n` bootstrap sampling of cells.
    """
    if verbose:
        print("fit TPT with source state: '" + source_state + "' and target state: '" + target_state + "'.")
    cluster_ident = {state: [int(i) for i in item] for state, item in cluster_ident.items()}
    flux = fit_tpt(transition_matrix, cluster_ident[source_state], cluster_ident[target_state], verbose = False)
    # netflux matrices of a 'null' background
    if verbose:
        print("determine transitions expected under a 'null' background.")
    random_flux = getReshuffledFlux(transition_matrix, cluster_ident,
                                    source_state, target_state, n = random_n)
    total_gross_flux_randomised = [m['tpt'].gross_flux.sum() for m in random_flux]
    (tpt_sets, tpt_coarse) = flux.coarse_grain([np.array(item) for c, item in cluster_ident.items()])
    # get the coarse-grained gross-flux and net-flux matrices
    gross_flux = tpt_coarse.gross_flux
    # normalised flux (i.e. % of flux over a given edge of the graph)
    gross_flux = gross_flux * 100.0/gross_flux.sum()

    # pathways
    (paths, pathfluxes) = tpt_coarse.pathways()
    path_dist = []
    # make sure the paths are named as in the user-defined coarse-graining scheme
    # this is the order in the paths printed
    cluster_order = [source_state] + [i for i in cluster_ident.keys() if i not in [source_state, target_state]] + [target_state]
    cluster_order = {n: j for n, j in enumerate(cluster_order)}
    for n, path in enumerate(paths):
        path = [cluster_order[j] for j in path]
        path_dist.append({'path': path, 'p_path': 100.0 * pathfluxes[n] / tpt_coarse.total_flux})
    
    # determine 'significance' of each transition (one-sided test,
    # i.e. number of times the observed gross-flux is larger than the randomised flux
    comparison = np.zeros(gross_flux.shape)
    random_gross_flux = []
    for m in random_flux:
        if type(m['gross_flux']) is np.ndarray:
            random_gross_flux.append(m['gross_flux'])
        else:
            random_gross_flux.append(m['gross_flux'].todense())
    random_gross_flux = np.array(random_gross_flux, dtype = 'float64')
    mean = np.mean(random_gross_flux[:, :, :])
    sd = np.std(random_gross_flux[:, :, :])
    for i in range(gross_flux.shape[0]):
        for j in range(gross_flux.shape[1]):
            comparison[i, j] = 1.0 - norm(loc = mean, scale = sd).cdf( gross_flux[i, j] )
    #for m in random_flux:
    #    comparison = np.add(comparison, np.greater_equal(gross_flux.todense(), m['gross_flux'].todense()))
    #comparison /= len(random_flux)
    
    # get bootstrap sampling of stationary distribution
    bootstrap_stationary = [getStationaryDistribution( flux.stationary_distribution, cluster_ident, i) for i in range(int(random_n))]
    bootstrap_stationary = { state: [j[state] for j in bootstrap_stationary] for state in cluster_ident.keys() }
    if type(gross_flux) is not np.ndarray:
        gross_flux = gross_flux.todense()

    return {'coarse_grain_tpt': tpt_coarse, 'gross_flux': gross_flux, \
            'gross_flux_randomised': random_gross_flux, \
            'pathways': pd.DataFrame(path_dist), 'randomised_tpt': random_flux, \
            'significance': comparison, 
            'total_gross_flux': tpt_coarse.gross_flux.sum(), \
            'total_gross_flux_randomised': total_gross_flux_randomised, \
            'stationary_bootstrapping': bootstrap_stationary }
