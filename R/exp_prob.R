# implement functions in python 'sceb' package for
# empirical Bayes estimator of gene-gene covariance matrix & inactive probabilities
# (as a better alternative for directly using the observed gene counts)

# all based on https://github.com/martinjzhang/single_cell_eb/blob/master/sceb/scdd.py
# (accessed 8 December 2021)

randomiseRowWise <- function(mat)
{
  # given a matrix, for each row randomise the columns
  t(apply(mat, MARGIN = 1, function(x){
    x[sample(1:length(x), length(x))]
  }))
}

dd_size_factor <- function(Y, assay.name = "IGHC")
{
  if( 'Seurat' %in% class( Y ) ){
    Y <- get_count_matrix( Y, assay.name = assay.name )
  } else if( 'dgCMatrix' %in% class( Y ) ){
    if( is.null( rownames( Y ) ) ) rownames( Y ) <- sapply(1:nrow(Y), function(y) paste0('g', y) )
    Y <- t(Y)
  }
  Nrc = Matrix::colSums( Y )
  Nr = mean( Nrc, na.rm = TRUE )
  size_factor = Nrc / Nr
  return( size_factor )
}

assign_row_weight <- function(X, row_weight)
{
  X <- Matrix::t(X) * row_weight
  return( Matrix::t(X) )
}

get_count_matrix <- function(SeuratObj, assay.name = "IGHC")
{
  return( SeuratObj@assays[[assay.name]]@counts )
}

get_JIC <- function(count_matrix)
{
  # modify the count matrix so that for each C gene we have counts for:
  # (1) IC (2) JC (3) IC + JC
  count_matrix <- count_matrix[ grepl( "-[IJ]C$", rownames( count_matrix ) ),  ]
  isotypes <- unique( sapply( rownames( count_matrix ), function(x) unlist(strsplit(x, split = "-"))[1] ) )
  out <- lapply( isotypes, function( isotype ){
    Matrix::colSums( count_matrix[ grepl( isotype, rownames( count_matrix ) ), ] )
  } )
  out <- do.call( "rbind", out )
  rownames( out ) <- paste0(isotypes, "-JIC")
  colnames( out ) <- colnames( count_matrix )
  out <- list( out, count_matrix[ grepl( "-JC$", rownames( count_matrix ) ),  ],
               count_matrix[ grepl( "-IC$", rownames( count_matrix ) ),  ] )
  return( do.call( "rbind", out ) )
}

dd_covariance <- function(SeuratObj, size_factor = NULL, assay.name = "IGHC", PC_prune = TRUE)
{
    # EB estimation of the covariance matrix and the Pearson correlation matrix.
    #Args:
    #    SeuratObj: Seurat Object with IGHC counts appended as a separate assay.
    #    size_factor (numeric vector, length Nc ): the cell size_factor.
    #    PC_prune (bool): if set the value to be zero for genes with small EB estimate variance (
    #        due to stability consideration)
    #Returns: a list of three items with names:
    #    'mean_dd' (numeric vector of length G): the mean gene expression level.
    #    'cov_dd' (matrix of dim [G, G]): the estimated covariance matrix.
    #    'PC_dd' (matrix of dim [G, G]): the estimated Pearson correlation matrix.
  X <- get_count_matrix( SeuratObj, assay.name = assay.name )
  G = dim( X )[1]; Nc <- dim( X )[2]
  gene_name <- rownames( X )

  # normalise by size_factor
  if( !is.null( size_factor ) ){
    row_weight <- 1 / size_factor
    X <- assign_row_weight( X, row_weight )
  }

  # Statistics after the first normalisation
  mean_dd <- rowMeans( X )
  M2_dd <- X %*% Matrix::t(X) / Nc

  # double normalisation to get what we need
  if( !is.null( size_factor ) ){
    X <- assign_row_weight( X, row_weight )
  }

  # Statistics after the second normalisation
  mean_doublenorm_dd <- rowMeans( X )

  # Bias correction
  Matrix::diag(M2_dd) <- Matrix::diag(M2_dd) - mean_doublenorm_dd

  # covariance matrix
  temp <- matrix( mean_doublenorm_dd )
  cov_dd <- M2_dd - temp %*% t(temp)
  diag_cov_dd <- Matrix::diag(cov_dd)
  index_bad <- rep(FALSE, G)
  index_bad[ which(diag_cov_dd <= 1e-1/9 )] <- TRUE
  index_bad[ which((diag_cov_dd / mean_dd) < 0.1) ] <- TRUE
  diag_cov_dd[ which(diag_cov_dd < 1e-12) ] <- 1e-12
  Matrix::diag(cov_dd) <- diag_cov_dd

  # Pearson correlation
  std_dd <- diag_cov_dd ^ 0.5
  temp <- matrix( std_dd )
  PC_dd <- cov_dd / ( std_dd %*% t(std_dd) )
  PC_dd[ PC_dd < -1 ] <- -1
  PC_dd[ PC_dd > 1 ] <- 1

  if( PC_prune ){
    PC_dd[, index_bad] <- 0
    PC_dd[index_bad, ] <- 0
  }

  return( list("mean_dd" = mean_dd, "cov_dd" = cov_dd, "PC_dd" = PC_dd) )
}

ml_covariance <- function(SeuratObj, size_factor = NULL, assay.name = "IGHC", PC_prune = TRUE)
{
  #  Plug-in estimation of the covariance matrix and the Pearson correlation matrix.
  #Args:
  #    SeuratObj: Seurat Object with IGHC counts appended as a separate assay.
  #    size_factor (numeric vector, length Nc ): the cell size_factor.
  #    PC_prune (bool): if set the value to be zero for genes with small EB estimate variance (
  #        due to stability consideration)
  #Returns: a list of three items with names:
  #    'mean_ml' (numeric vector of length G): the mean gene expression level.
  #    'cov_ml' (matrix of dim [G, G]): the estimated covariance matrix.
  #    'PC_ml' (matrix of dim [G, G]): the estimated Pearson correlation matrix.
  X <- get_count_matrix( SeuratObj, assay.name = assay.name )
  G = dim( X )[1]; Nc <- dim( X )[2]
  Nr <- sum( X ) / Nc
  gene_name <- rownames( X )

  # normalise by size_factor
  if( !is.null( size_factor ) ){
    row_weight <- 1 / size_factor
    X <- assign_row_weight( X, row_weight )
  }

  # mean
  mean_ml <- rowMeans( X )

  # 2nd moment
  M2_ml <- X %*% Matrix::t(X) / Nc

  # covariance matrix
  temp <- matrix( mean_ml )
  cov_ml <- M2_ml - temp %*% t(temp)

  # Pearson correlation
  std_ml <- Matrix::diag(cov_ml)
  std_ml[ which(std_ml < 0) ] <- 0
  std_ml <- std_ml ^ 0.5
  temp <- matrix( std_ml )
  PC_ml <- cov_ml / temp %*% t(temp)
  PC_ml[ PC_ml < -1 ] <- -1
  PC_ml[ PC_ml > 1 ] <- 1

  return( list("mean_ml" = mean_ml, "cov_ml" = cov_ml, "PC_ml" = PC_ml) )

}

smooth_zero_estimator <- function(L, t = 1, n = 500, require_param = FALSE, restrict_t = TRUE)
{
  L <- ceiling( L )
  if( restrict_t ){
    t <- min( t, 5) # for robustness consideration
  }
  v_L <- seq( from = 0, to = round(L) )
  t_ <- t - 1
  w <- (- t_) ** v_L
  if( t > 2 ){
    k <- ceiling( 0.5 * log2( n * t_ ** 2 / (t_ - 1)) )
    w <- w * ( 1 - stats::pbinom(v_L - 1, k, 1 / t))
  }
  if( require_param ){
    return( list('w' = w, 't' = t, 'k' = k, 'q' = 1 / (t + 1)) )
  } else {
    return( w )
  }
}

bincount <- function(x, n_gene, weights = 1)
{
  # implementation to do the equivalent of numpy.bincount
  # JN (18-01-2022): the original is very slow in MacOS,
  # use plyr::ddply to speed up summing the weights
  lvls <- seq(from = 1, n_gene)
  if( weights != 1 && length( x ) != length( weights ) )
    stop("Lengths of 'x' and 'weights' in bincount do not match.")
  #sapply( 1:length(lvls), function(n){
  #  sum(as.numeric(x == lvls[n]) * weights)
  #})
  if( length(weights) == 1 && weights == 1 ) weights <- rep(1, length(x))
  dt <- data.frame(vec = x, w = weights)
  dt$vec <- factor(dt$vec, levels = lvls)
  rez <- plyr::ddply(dt, .variables = "vec", function(x) sum(x[, "w"]))#plyr::summarise, wg_sum = sum(w), .drop = FALSE)
  o <- as.vector(rez[, 2])
  names(o) <- rez[, 1]
  o <- o[lvls]
  o[is.na(o)] <- 0
  unname(o)
}

dd_inactive_prob <- function(Y, genes = NULL, relative_depth = 1, size_factor = NULL, assay.name = 'IGHC')
{
  sub_dd_inactive_prob <- function(Y_sub, t_sub, Nc_sub, n = 500, n_genes = G)
  {
    # Y_sub is a dgCMatrix
    # A is the vector of non-zero entries of Y_sub
    A = Y_sub@x
    # column indices of each observation in Y_sub@x (= A)
    JA = Y_sub@i + 1# findInterval(seq(Y_sub@x) - 1, Y_sub@p[-1]) + 1
    if( t_sub < 1 ){
      L <- log(0.005) / log(1 - t_sub)
    } else{
      L <- 10
    }
    w <- smooth_zero_estimator( L, t = t_sub, n = n, require_param = FALSE )
    p0_ml <- rep(0, n_genes); p0_dd <- rep(0, n_genes)
    A_dd <- rep(0, length(A))
    for(i in 1:length(w)){
      A_dd[which(A == (i - 1))] <- w[i]
    }
    # use bincount above to substitute np.bincount in python.
    p0_ml = 1 - bincount(JA, n_gene = n_genes) / Nc_sub
    p0_dd = p0_ml + bincount(JA, n_gene = n_genes, weights = A_dd) / Nc_sub
    return( list('p0_ml' = p0_ml, 'p0_dd' = p0_dd) )
  }
  if( 'Seurat' %in% class( Y ) ){
    Y <- get_count_matrix( Y, assay.name = assay.name )
  } else if( 'dgCMatrix' %in% class( Y ) ){
    if( is.null( rownames( Y ) ) ) rownames( Y ) <- sapply(1:nrow(Y), function(y) paste0('g', y) )
  }
  if( is.null( genes ) ){
    genes <- rownames( Y )
  } else {
    # check whether all items listed in 'genes' is in rownames( Y )
    not_in <- genes[ which(! genes %in% rownames( Y ))]
    if( length( not_in ) > 0 ){
      stop( paste0( "A total of ", length( not_in ), " genes are not found in the indicated assay. Make sure you have supplied the correct SeuratObj, assay.name and genes." ) )
    }
  }
  Y <- Y[ genes, which(Matrix::colSums(Y) > 0) ]
  G = dim( Y )[1]; Nc <- dim( Y )[2]
  # Y <- Matrix::t( Y )
  n <- min(500, Nc)
  Y <- as.matrix ( Y )# just make the subsetting easier
  if( is.null(size_factor) ){
    p0 <- sub_dd_inactive_prob( Y, 1 / relative_depth, Nc, n = n, n_genes = G)
    p0_ml <- p0[['p0_ml']]; p0_dd <- p0[['p0_dd']]
  } else {
    size_factor <- size_factor[ colnames(Y) ]
    tc <- 1 / relative_depth / size_factor
    amp <- max(20, 1 / stats::quantile(tc, probs = 0.001))
    tc <- round( tc * amp ) / amp
    tc[ which(tc < 1 / amp)] <- 1 / amp
    p0_ml <- rep(0, G); p0_dd <- rep(0, G)
    for(tc_ in sort(unique(tc))){
      Nc_sub <- sum(tc == tc_)
      if( Nc_sub > 0 ){
        Y_sub <- Matrix::Matrix( Y[, which(tc == tc_)], sparse = TRUE )
        p0_sub <- sub_dd_inactive_prob( Y_sub, tc_, Nc_sub, n = n )
        p0_ml_sub <- p0_sub[['p0_ml']]; p0_dd_sub <- p0_sub[['p0_dd']]
        p0_ml <- p0_ml + p0_ml_sub * Nc_sub / Nc
        p0_dd <- p0_dd + p0_dd_sub * Nc_sub / Nc
      }
    }
  }
  p0_ml[ p0_ml < 0 ] <- 0
  p0_dd[ p0_dd < 0 ] <- 0
  return( list( 'p0_ml' = p0_ml, 'p0_dd' = p0_dd ) )
}

dd_pairwise_inactive_prob <- function(Y, genes = NULL, relative_depth = 1, size_factor = NULL, assay.name = 'IGHC')
{
  sub_dd_pairwise_inactive_prob <- function(Y_sub, t_sub, Nc_sub, n = 500, n_genes = G)
  {
    # Y_sub is a dgCMatrix
    # row indices of each observation in Y_sub@x (= A)
    IA = Y_sub@i + 1
    # column indices of each observation in Y_sub@x (= A)
    JA = Y_sub@j + 1# findInterval(seq(Y_sub@x) - 1, Y_sub@p[-1]) + 1
    A = Y_sub@x
    if( t_sub < 1 ){
      L <- log(0.005) / log(1 - t_sub)
    } else{
      L <- 10
    }
    #print(paste0("L: ", L))
    w <- smooth_zero_estimator( L, t = t_sub, n = n, require_param = FALSE )
    #print(paste0("w: ", paste(head(w), collapse = " ")))
    # Maintain weights for the dd estimator. The weights of ml is always 1, no need of maintenance
    A_dd <- rep(0, length(A))
    for(i in 1:length(w)){
      A_dd[which(A == (i - 1))] <- w[i]
    }
    #print(head(A_dd))
    zero_matrix_dd <- matrix(0, nrow = n_genes, ncol = n_genes)
    zero_matrix_ml <- matrix(0, nrow = n_genes, ncol = n_genes)

    # use bincount above to substitute np.bincount in python.
    #print(paste0("JA: ", paste(head(JA), collapse = " ")))
    temp_J_list <- bincount(JA, n_gene = n_genes)
    temp_w_list <- bincount(JA, n_gene = n_genes, weights = A_dd)
    for(i_gene in 1:length( temp_J_list )){
      #i_gene <- i_gene - 1
      zero_matrix_ml[i_gene, ] <- zero_matrix_ml[i_gene, ] + temp_J_list[i_gene]
      zero_matrix_ml[, i_gene] <- zero_matrix_ml[, i_gene] + temp_J_list[i_gene]
      zero_matrix_dd[i_gene, ] <- zero_matrix_dd[i_gene, ] + temp_w_list[i_gene]
      zero_matrix_dd[, i_gene] <- zero_matrix_dd[, i_gene] + temp_w_list[i_gene]
    }

    # update the intersection part of the ml matrix
    temp_ml <- Matrix::sparseMatrix(x = rep(1, length(A)), i = IA, j = JA,
                                    dims = c(Nc_sub, G))
    zero_matrix_ml <- zero_matrix_ml - as.matrix( Matrix::t( temp_ml ) %*% temp_ml )

    # update the intersection part of the dd matrix
    temp_dd <- Matrix::sparseMatrix(x = A_dd, i = IA, j = JA, dims = c(Nc_sub, G))
    zero_matrix_dd <- zero_matrix_dd + as.matrix( Matrix::t( temp_dd ) %*% temp_dd )
    temp <- as.matrix( Matrix::t( temp_ml ) %*% temp_dd )
    zero_matrix_dd <- zero_matrix_dd - temp - t( temp )

    zero_matrix_ml <- 1 - zero_matrix_ml / Nc_sub
    zero_matrix_dd <- zero_matrix_dd / Nc_sub + zero_matrix_ml
    #print(zero_matrix_ml[1:3, 1:3])
    #print(zero_matrix_dd[1:3, 1:3])
    return( list('p0_ml' = zero_matrix_ml, 'p0_dd' = zero_matrix_dd) )
  }
  # first calculate diagonal elements and set aside
  p0 <- dd_inactive_prob(Y, relative_depth = relative_depth, size_factor = size_factor,
                         assay.name = assay.name, genes = genes)
  if( 'Seurat' %in% class( Y ) ){
    Y <- get_count_matrix( Y, assay.name = assay.name )
  } else if( 'dgCMatrix' %in% class( Y ) ){
    if( is.null( rownames( Y ) ) ) rownames( Y ) <- sapply(1:nrow(Y), function(y) paste0('g', y) )
  }
  if( is.null( genes ) ){
    genes <- rownames( Y )
  } else {
    # check whether all items listed in 'genes' is in rownames( Y )
    not_in <- genes[ which(! genes %in% rownames( Y ))]
    if( length( not_in ) > 0 ){
      stop( paste0( "A total of ", length( not_in ), " genes are not found in the indicated assay. Make sure you have supplied the correct SeuratObj, assay.name and genes." ) )
    };
  }
  Y <- Y[ genes,  which(Matrix::colSums(Y) > 0)]
  G = dim( Y )[1]; Nc <- dim( Y )[2]
  # Y <- Matrix::t( Y )
  n <- min(500, Nc)
  Y <- as.matrix ( Y )# just make the subsetting easier
  if( is.null(size_factor) ){
    zero_matrix <- sub_dd_pairwise_inactive_prob( Y, 1 / relative_depth, Nc, n = n, n_genes = G)
    zero_matrix_ml <- zero_matrix[['p0_ml']]; zero_matrix_dd <- zero_matrix[['p0_dd']]
  } else {
    size_factor <- size_factor[ colnames(Y) ]
    tc <- 1 / relative_depth / size_factor
    amp <- max(20, 1 / stats::quantile(tc, probs = 0.001))
    tc <- round( tc * amp ) / amp
    tc[ which(tc < 1 / amp)] <- 1 / amp
    zero_matrix_ml <- matrix(0, nrow = G, ncol = G)
    zero_matrix_dd <- matrix(0, nrow = G, ncol = G)
    # set a progress bar
    #pb <- utils::txtProgressBar( min = 0, max = length( unique( tc ) ),
    #                             style = 3, file = stderr() )
    i <- 1
    for(tc_ in sort(unique(tc))){
      Nc_sub <- sum(tc == tc_)
      #print(paste0("tc_: ", tc_))
      if( Nc_sub > 0 ){
        Y_sub <- Matrix::Matrix( Y[, which(tc == tc_)], sparse = TRUE)
        Y_sub <- Matrix::t( Y_sub )
        Y_sub <- methods::as(Y_sub, "dgTMatrix")
        zero_matrix_sub <- sub_dd_pairwise_inactive_prob( Y_sub, tc_, Nc_sub, n = n )
        zero_matrix_sub_ml <- zero_matrix_sub[['p0_ml']]
        zero_matrix_sub_dd <- zero_matrix_sub[['p0_dd']]
        zero_matrix_ml <- zero_matrix_ml + zero_matrix_sub_ml * Nc_sub / Nc
        zero_matrix_dd <- zero_matrix_dd + zero_matrix_sub_dd * Nc_sub / Nc
      }
      #utils::setTxtProgressBar( pb = pb, value = i)
      i <- i + 1
    }
  }

  # fill in diagonals
  diag( zero_matrix_ml ) <- p0[['p0_ml']]
  diag( zero_matrix_dd ) <- p0[['p0_dd']]

  zero_matrix_ml[ zero_matrix_ml < 0 ] <- 0
  zero_matrix_dd[ zero_matrix_dd < 0 ] <- 0
  colnames(zero_matrix_ml) <- genes
  colnames(zero_matrix_dd) <- genes
  rownames(zero_matrix_ml) <- genes
  rownames(zero_matrix_dd) <- genes
  return( list( 'p0_ml' = zero_matrix_ml, 'p0_dd' = zero_matrix_dd ) )
}

pairwise_zero_to_coexpression <- function(mat)
{
  # convert a pairwise-zero-expression matrix to a pairwise co-expression matrix by:
  # p'(1,2) = 1 - p(1,2) - p1 * (1 - p2) - p2 * (1 - p1)
  # where p1, p2 are probability of zero expression of genes 1 & 2,
  #       p(1,2) is the pairwise zero probabiility of genes 1 & 2
  #       p'(1,2) is the pairwise coexpression probability of genes 1 & 2
  dn <- dimnames( mat )
  o <- sapply(1:nrow(mat), function(i){
    sapply(1:ncol(mat), function(j){
      if(i != j){
        a <- c(mat[i, i], mat[j, j])
        o <- 1 - mat[i, j] - (a[1] * (1 - a[2]) + a[2] * (1 - a[1]))
        if( o < 0 ) o <- 0
        return(o)
      } else return(1 - mat[i, j])
    })
  })
  dimnames( o ) <- dn
  return( o )
}

get_coex_by_cluster <- function( SeuratObj, assay.name = 'IGHC', group.by = 'seurat_clusters' )
{
  size_factor <- dd_size_factor( SeuratObj, assay.name = assay.name )
  count_matrix <- get_count_matrix( SeuratObj, assay.name = assay.name )
  count_matrix <- get_JIC( count_matrix )
  genes <- rownames( count_matrix )
  cluster_info <- Seurat::FetchData( SeuratObj, group.by )
  if( is.factor( cluster_info[, group.by] ) ){
    cluster_list <- levels( cluster_info[, group.by] )
  } else {
    cluster_list <- sort( unique( cluster_info[, group.by] ) )
  }
  o <- lapply( cluster_list, function( cluster ){
    y <- count_matrix[, which( colnames( count_matrix ) %in% rownames( cluster_info )[ which( cluster_info[, group.by] == cluster) ] )]
    message('\nCalculating cluster ', cluster, ' ...')
    csr_history <- dd_pairwise_inactive_prob( y, genes = genes[grepl("-JC$|-JIC$", genes)], size_factor = size_factor[ colnames(y) ])$p0_dd
    csr_potential <- dd_pairwise_inactive_prob( y, genes = genes[grepl("-JC$|-IC$", genes)], size_factor = size_factor[ colnames(y) ])$p0_dd
    csr_history <- pairwise_zero_to_coexpression( csr_history )
    csr_potential <- pairwise_zero_to_coexpression( csr_potential )
    csr_history <- reshape2::melt( csr_history )
    csr_potential <- reshape2::melt( csr_potential )
    colnames( csr_history ) <- c("JC", "int", "proba") # 'int' for intermediates (of productive transcripts)
    colnames( csr_potential ) <- c("JC", "int", "proba")
    csr_history$partition <- 'CSR history'
    csr_potential$partition <- 'CSR potential'
    csr_history <- csr_history[ which(grepl('-JC$', csr_history$JC) &
                                        grepl('IC$', csr_history$int)), ]
    csr_potential <- csr_potential[ which(grepl('-JC$', csr_potential$JC) &
                                            grepl('IC$', csr_history$int)), ]
    ic_jc_co <- rbind( csr_history, csr_potential )
    ic_jc_co$int <- gsub("-[JI]*C$", "", gsub("^IGH", "", ic_jc_co$int))
    ic_jc_co$JC <- gsub("-[JI]*C$", "", gsub("^IGH", "", ic_jc_co$JC))
    ic_jc_co$int <- factor(ic_jc_co$int, levels = c("M", "D", "G3", "G1",
                                                    "A1", "G2", "G4", "E", "A2"),
                          ordered = TRUE)
    ic_jc_co$JC <- factor(ic_jc_co$JC, levels = c("M", "D", "G3", "G1",
                                                  "A1", "G2", "G4", "E", "A2"),
                          ordered = TRUE)
    ic_jc_co$int_n <- as.numeric(ic_jc_co$int)
    ic_jc_co$JC_n <- as.numeric(ic_jc_co$JC)
    ic_jc_co$partition <- apply(ic_jc_co[, c("int_n", "JC_n", "partition")], MARGIN = 1, function(x){
      if(x[1] == x[2]) return('diag')
      ifelse( as.numeric(x[1]) < as.numeric(x[2]), 'CSR history', 'CSR potential')
    })
    ic_jc_co <- ic_jc_co[ ic_jc_co$partition != 'diag', ]
    ic_jc_co <- ic_jc_co[ , c('JC', 'int', 'partition', 'proba') ]
    ic_jc_co[ , group.by ] <- cluster
    return( ic_jc_co )
  } )
  return( do.call( "rbind", o ) )
}
