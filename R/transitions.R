#' read loom matrices
#'
#' @description
#' reproduce read.loom.matrices from velocyto.R to avoid problems in installing dependencies
#'
#' @param file input loom file
#' @param engine package to read h5 file type (only 'hdf5r' supported now)
#'
#' @importFrom hdf5r H5File list.datasets
#' @importFrom methods as
#' @export read_loom_matrices
read_loom_matrices <- function (file, engine = "hdf5r")
{
  if (engine == "hdf5r") {
    # cat("reading loom file via hdf5r...\n")
    f <- hdf5r::H5File$new(file, mode = "r")
    cells <- f[["col_attrs/CellID"]][]
    genes <- f[["row_attrs/Gene"]][]
    dl <- c(spliced = "layers/spliced", unspliced = "layers/unspliced",
            ambiguous = "layers/ambiguous")
    if ("layers/spanning" %in% hdf5r::list.datasets(f)) {
      dl <- c(dl, c(spanning = "layers/spanning"))
    }
    dlist <- lapply(dl, function(path) {
      m <- as(t(f[[path]][, ]), "dgCMatrix")
      rownames(m) <- genes
      colnames(m) <- cells
      return(m)
    })
    f$close_all()
    return(dlist)
  }
  else {
    warning("Unknown engine. Use hdf5r to import loom file.")
    return(list())
  }
}

#' get somatic hypermutation level
#'
#' @description
#' \code{getSHM} calculates somatic hypermutation (i.e. 1 - (percentage identity to germline VH gene)).
#'
#' @details
#' \code{getSHM} considers the \code{v_identity} (i.e. % identity to germline V gene) for the VH sequence
#' mapped to each cell, and calculate 1 - \code{v_identity} as the somatic hypermutation (SHM)
#' level for the cell. It finds the v_identity information from the meta.data slot of the Seurat object.
#' For cells without a mapped V sequence, it will impute SHM as 0.
#'
#' @param SeuratObj Seurat Object
#' @param v_identity_anno_name column name in the Seurat meta.data slot which holds the v_identity information to be considered.
#' @param shm_column_to_add name of column to be added to the Seurat meta.data slot
#' which holds the SHM frequency calculated in this function
#'
#' @return The same Seurat object as given by SeuratObj, except that a new column with
#' name given by \code{shm_column_to_add} is appended to the SeuratObj meta.data slot to reflects
#' the calculated SHM frequency.
#'
#' @importFrom Seurat AddMetaData
#' @export getSHM
getSHM <- function(SeuratObj, v_identity_anno_name,
                   shm_column_to_add = "shm")
{
  v_identity <- SeuratObj[[v_identity_anno_name]]
  shm <- 1 - v_identity[, v_identity_anno_name]
  shm <- (shm - min(shm, na.rm = TRUE)) / abs(diff(range(shm, na.rm = TRUE)))
  shm[is.na(shm)] <- 0
  v_identity[, shm_column_to_add] <- shm
  Seurat::AddMetaData(SeuratObj, v_identity[, shm_column_to_add, drop = FALSE])
}

#' score cells by their Class Switch Recombination (CSR) status
#'
#' @description
#' \code{getCSRpotential} scores each cell by their status in temrs of class switch recombination (CSR), by considering the mapped sterile and productive IgH transcripts.
#'
#' @details
#' \code{getCSRpotential} calculates a "CSR potential" score which ranks the cells in the given Seurat Obj by their status in the CSR process. The default is to calculate this by estimating the contribution (weight) of a 'naive' isotype signature for each cell, given its sterile/productive expression profile.
#' The CSR potential will be 1 - (Naive signature weight). This method is available for either human or mouse for which isotype signatures were trained on reference B cell atlas data. Alternatively, CSR potential can also be calculated empirically (by setting \code{reference_based = NULL}),
#' given by the Euclidean norm of (representative_p, total_s) (i.e.\eqn{ \sqrt{ \text{representative_p}^2 + \text{total_s}^2} } ), where
#' \itemize{
#'   \item{\code{representative_p}}{: the productive isotype for each cell (in human this will be 0 = IgM, 1 = IgG3 ... ), and}
#'   \item{\code{total_s}}{: amount of sterile IgH molecules for each cell.}
#' }
#' For \code{total_s}, the default is to use the scale.data slot which already normalises the IGHC counts by library size. If this doesn't exist the function will calculate this while regressing out the library size.
#' For \code{representative_p}, users can either use a specified column in the Seurat object meta.data which indicates the isotype of the cell, or, if not provided, used the productive reads counted using the productive/sterile quantification workflow implemented in this package (see the argument \code{mode} of this function).
#'
#' @param SeuratObj Seurat Object
#' @param ighc_count_assay_name name of assay in \code{SeuratObj} which holds the IgH productive/sterile transcript count data. (Default: "IGHC")
#' @param ighc_slot the slot in \code{slot(SeuratObj, "assays")[[ighc_count_assay_name]]} to be used to access productive/sterile transcript counts (Default: "scale_data")
#' @param knn_graph should the k-nearest neighbour graph calculated on the gene expression assay be used to impute the annotation of productive transcripts for cells where no such transcripts are found across all isotypes? If TRUE, majority voting on the direct neighbours of the cell in the kNN graph will be used to impute. Otherwise, the cell will be assume to express IgM productive transcript. Expects \code{TRUE} or \code{FALSE}, or a \code{igraph} object containing kNN graph (in which case this graph will be used for majority voting imputation). (Default: TRUE)
#' @param reference_based indicate the species. The function will use a naive isotype signature (sterile/productive gene counts) trained on reference B cell atlas for the given species. For now either 'human' or 'mouse' are accepted. If \code{NULL}, the function calculates CSR potential by taking the Euclidean norm of (representative_p, total_s) (see Details).
#' @param vars.to.regress list of variables to be regressed out when scaling the sterile/productive count matrix, if \code{ighc_slot} is given as \code{scale.data} but it has not been populated. (Default: "nCount_RNA", i.e. per-cell library size)
#' @param mode (Only applicable if c_gene_anno_type is NULL.) Interpretation of the isotype expressed by the cell. Either "furthest" (i.e. the isotype furthest along the IGH locus with non-zero expression of productive transcript will be taken as the isotype representative of the cell) or "highest" (the isotype with highest expression). (Default: "furthest")
#' @param c_gene_anno_name If not NULL, this column from the Seurat Object meta.data will be used to indicate \code{representative_p} in calculaing the CSR potential score, in lieu of the productive transcript counts in the IGHC assay (Default: \code{NULL})
#' @param isotype_column_to_add name of column to be added to the SeuratObj meta.data to indicate the isotype of the cell. Used for subsequent grouping of cells in calculating transitions.
#'
#' @return Seurat object with these following columns added to the meta.data slot:
#' \itemize{
#'   \item{\code{representative_p}: }{an integer indicating the productive isotype for each cell (in human: 0 = IgM, 1 = IgG3 ... )}
#'   \item{\code{total_s}: }{amount of sterile IgH molecules for each cell, calculated from the given \code{ighc_slot} of the IGHC assay.}
#'   \item{\code{csr_pot}: }{CSR potential. Depending on the argument \code{reference_based} the method of calculation will be different (see Details).}
#'   \item{\code{isotype_column_to_add}: }{isotype labelled as M, G3, etc. (added only when c_gene_anno_name is FALSE and the ighc_count_assay_name Assay is used to calculate CSR potential.}
#' }
#'
#' @importFrom Seurat Assays ScaleData AddMetaData
#' @importFrom Matrix colSums
#' @importFrom igraph graph_from_adjacency_matrix neighbors
#' @importFrom methods slotNames slot
#' @importFrom nnls nnls
#' @importFrom stringr str_replace str_to_sentence str_to_upper
#' @export getCSRpotential
getCSRpotential <- function(SeuratObj, ighc_count_assay_name = "IGHC",
                            ighc_slot = "scale.data", knn_graph = TRUE,
                            reference_based = NULL,
                            vars.to.regress = c("nCount_RNA"),
                            mode = "furthest", c_gene_anno_name = NULL,
                            isotype_column_to_add = "isotype")
{
  if(! ighc_count_assay_name %in% Seurat::Assays(SeuratObj))
    stop(paste0("The assay '", ighc_count_assay_name, "' cannot be found in SeuatObj.") )
  if( ! ighc_slot %in% slotNames(SeuratObj@assays[[ighc_count_assay_name]]) )
    stop(paste0("The named slot '", ighc_slot, "' cannot be found in the '",
                ighc_count_assay_name, "' assay in SeuratObj."))
  if( ! mode %in% c("furthest", "highest") )
    stop(paste0("'mode' must be either 'furthest' or 'highest'."))
  if( !is.null( reference_based) ){
    if( ! reference_based %in% c("human", "mouse") )
      stop("'reference_based' must be either 'human' or 'mouse'. If your data come from another species, or you wish to calculate CSR potential empirically, set reference_based = NULL. ")
  }
  if( !is.logical( knn_graph ) ){
    if( inherits( knn_graph, 'igraph') ){
      stop("'knn_graph' must either be TRUE or FALSE, or an igraph object of the k-nearest neighbour graph.")
    }
  } else {
    if( knn_graph ){
      default_assay <- Seurat::DefaultAssay( SeuratObj )
      if( ! paste0(default_assay, "_nn") %in% names( slot(SeuratObj, "graphs") ) ){
        stop("no k-nearest neighbour graph found. Please run Seurat::FindNeighbors first.")
      }
      knn_graph <- igraph::graph_from_adjacency_matrix( slot(SeuratObj, "graphs")[[ paste0(default_assay, "_nn") ]],
                                                        diag = FALSE, mode = "undirected" )
    }
  }

  if( ighc_slot == "scale.data" &&
      all(dim(SeuratObj@assays[[ighc_count_assay_name]]@scale.data) == 0)){
    # populate the scale.data matrix by running Seurat::ScaleData
    SeuratObj <- Seurat::ScaleData(SeuratObj, vars.to.regress = vars.to.regress,
                                   assay = ighc_count_assay_name)
  }
  ighc_counts <- slot(SeuratObj@assays[[ighc_count_assay_name]], ighc_slot)
  raw_counts <- slot(SeuratObj@assays[[ighc_count_assay_name]], "counts")
  ic_genes <- rownames(ighc_counts)[grepl("S$", rownames(ighc_counts))]
  total_ic <- Matrix::colSums(ighc_counts[ic_genes, ])
  c_genes <- rownames(ighc_counts)[grepl("-C$", rownames(ighc_counts))]
  c_genes <- stringr::str_replace(stringr::str_replace(c_genes, "-C$", ""),
                                  "^IGH|^Igh", "")
  c_genes <- stringr::str_to_sentence(c_genes)

  if( is.null( c_gene_anno_name ) ){
    if( is.null( isotype_column_to_add ) )
      stop("'isotype_column_to_add' needs to be given if 'c_gene_anno_name' is NULL.")
    jc_genes <- rownames(ighc_counts)[grepl("P$", rownames(ighc_counts))]
    jc_counts <- raw_counts[jc_genes, ] # use raw counts for productive reads
    furthest_jc <- sapply(colnames(jc_counts), function(y) {
      x <- jc_counts[, y]
      if( mode == "furthest" ) {
	      pos <- which(x > 0) # positive entries
        if(length(pos) == 0) { # no productive IgH detected
	        if( inherits(knn_graph, 'igraph') ){
	          # majority voting using neighbours in the kNN graph where this information does exist
	          neighbours <- names(igraph::neighbors(knn_graph, y))
	          if( length(neighbours) == 0 ) return( 0 ) # assume IgM
	          neighbours_jc <- apply(jc_counts[, neighbours, drop = FALSE], 2, function(xx){
	            poss <- which(xx > 0)
	            if( length(poss) > 0 ) return(max(poss) - 1) else return(NA)
            })
	          o <- as.numeric( names(which.max(table(neighbours_jc))) )
	          if( length( o ) == 0 ) return(0) # assume IgM
	          else return(o)
	        } else return(0) # assume IgM
	      } else return(max(pos) - 1)
      } else if ( mode == "highest" ){
        pos <- which(x > 0) # positive entries
        if(length(pos) == 0) { # no productive IgH detected
          if( inherits(knn_graph, 'igraph') ){
            # majority voting using neighbours in the SNN graph where this information does exist
            neighbours <- names(igraph::neighbors(knn_graph, y))
	    if( length(neighbours) == 0 ) return( 0 ) # assume IgM
            neighbours_jc <- apply(jc_counts[, neighbours, drop = FALSE], 2, function(xx){
              poss <- which(xx > 0)
              if( length(poss) > 0 ) return(max(which(xx == max(xx, na.rm = TRUE))) - 1) else return(NA)
            })
            o <- as.numeric( names(which.max(table(neighbours_jc))) )
	    if( length( o ) == 0 ) return(0) else return( o )
          } else return(0) # assume IgM
        } else return(max(which(x == max(x, na.rm = TRUE))) - 1)
      }
    })
    # now change these numbers to string and add to SeuratObj meta.data
    c_gene_anno <- factor(furthest_jc, levels = (1:length(c_genes)) - 1,
                          labels = c_genes)
    names(c_gene_anno) <- colnames(ighc_counts)
    SeuratObj <- Seurat::AddMetaData(SeuratObj, c_gene_anno,
                                     col.name = isotype_column_to_add)
  } else {
    if( ! c_gene_anno_name %in% colnames(SeuratObj@meta.data))
      stop("The column name given in 'c_gene_anno_name' is not found in SeuatObj@meta.data." )
    jc_anno <- SeuratObj@meta.data[, c_gene_anno_name]
    # create a numeric order of jc_anno. Trouble is it can take different format
    # ('IGHM', 'IgM', 'M', 'Ighm' etc.)
    # The implementation: first remove all 'IGH', 'Ig', 'Igh' etc. first, and then
    # convert the remaining substring into upper case.
    # use the order in ighc_count c_genes to generate a factor and from there get the numerical order
    jc_anno <- stringr::str_replace(jc_anno, "^IGH|^Igh", "")
    jc_anno <- stringr::str_replace(jc_anno, "^Ig", "")
    jc_anno <- stringr::str_to_upper(jc_anno)
    c_genes <- stringr::str_replace(stringr::str_replace(ic_genes, "-S$", ""),
                                    "^IGH|^Igh", "")
    jc_anno <- factor(jc_anno, levels = stringr::str_to_upper(c_genes))
    furthest_jc <- as.numeric(jc_anno) - 1
    names(furthest_jc) <- Seurat::Cells(SeuratObj)
    if( inherits(knn_graph, 'igraph') ){
      furthest_jc <- sapply(names(furthest_jc), function(x){
        if( !is.na(furthest_jc[x]) ) return( furthest_jc[x] )
        else{
          # majority voting using neighbours in the SNN graph where this information does exist
          neighbours <- names(igraph::neighbors(knn_graph, x))
	  if( length(neighbours) == 0 ) return(0) # assume IgM
          neighbours_jc <- furthest_jc[neighbours]
          o <- as.numeric( names(which.max(table(neighbours_jc))) )
          if( length( o ) == 0 ) return(0) else return(o)
	}
      })
    } else {
      furthest_jc[is.na(furthest_jc)] <- 0
    }
    c_gene_anno <- factor(furthest_jc, levels = (1:length(c_genes)) - 1,
                          labels = stringr::str_to_sentence(c_genes))
    names(c_gene_anno) <- colnames(ighc_counts)
    SeuratObj <- Seurat::AddMetaData(SeuratObj, c_gene_anno,
                                     col.name = isotype_column_to_add)
  }
  if( !is.null(reference_based) && reference_based %in% c("human", "mouse") ){
    # use a reference 'signature' matrix and apply non-negative least-square
    # to infer the coefficients of each signature given the expression profile
    # of Igh sterile/productive transcripts
    # CSR potential = 1 - coef of the naive signature
    assign("nmf_signatures", get(paste0(reference_based, "_nmf")))
    naive_signature <- sapply(Cells(SeuratObj),
                              function(n) nnls::nnls(nmf_signatures,
                                                     ighc_counts[rownames(nmf_signatures), n])$x[1])
    # scale into range [0, 1] and invert
    naive_signature <- (naive_signature - min(naive_signature)) / abs(diff(range(naive_signature)))
    csr_pot <- 1 - naive_signature
  } else {
    # use the empirical method = sqrt( furthest_productive ^ 2 + total_sterile ^ 2)
    csr_pot <- sqrt( furthest_jc^2 + total_ic^2)
    csr_pot <- (csr_pot - min(csr_pot, na.rm = TRUE)) / abs(diff(range(csr_pot, na.rm = TRUE)))
    furthest_jc[is.na(furthest_jc)] <- 0
    total_ic[is.na(total_ic)] <- 0
    csr_pot[is.na(csr_pot)] <- 0
  }
  o <- data.frame(representative_p = furthest_jc, total_s = total_ic,
                  csr_pot = csr_pot)
  rownames(o) <- colnames(ighc_counts)
  SeuratObj <- Seurat::AddMetaData(SeuratObj, o)
  return(SeuratObj)
}

#' wrapper function to convert Seurat Object to a AnnData .h5ad file
#'
#' @description
#' \code{convertSeuratToH5ad} is a wrapper function to convert a given Seurat Object into an AnnData object (for use in python with e.g. scanpy) and write out into a .h5ad file.
#'
#' @details
#' \code{convertSeuratToH5ad} simply wraps around the R \code{SeuratDisk} package to perform the stated conversion. Included her for ease of use for the user.
#' Each assay in the Seurat Object is written into separate .h5ad files.
#'
#' @param SeuratObj Seurat Object
#' @param assays A vector of assay names from SeuratObj to be exported
#' @param h5ad_filename Filename of the output .h5ad file. Note that the final output filenames will have the assay names appended to this (see examples).
#' @param conda_env character, if not \code{NULL} this named conda environment is used.
#' (Default: 'scicsr'). If \code{NULL}, no conda environment will be used, the program assumes the python packages \code{scanpy} and \code{scvelo} are installed in the local python)
#'
#' @return A vector of .h5ad filenames which are outputted. Each file correspond to one Seurat assay, as indicated in the suffix inside the filename (see examples).
#'
#' @import reticulate
#' @importFrom Seurat Assays
#' @importFrom SeuratDisk SaveH5Seurat Convert
#' @export convertSeuratToH5ad
convertSeuratToH5ad <- function(SeuratObj, assays, h5ad_filename,
                                conda_env = 'scicsr'){
  # first convert all metadata columns which are factors into characters
  # this is to preserve them in the anndata object (if they are factors only the levels
  # will be retained, not the labels)
  for(col in colnames(SeuratObj@meta.data)){
    if(inherits(SeuratObj@meta.data[, col], "factor")){
      SeuratObj@meta.data[, col] <- as.character( SeuratObj@meta.data[, col] )
    }
  }
  if( ! grepl(".h5ad$", h5ad_filename) )
    stop("The filename given in 'h5ad_filename' should end with the extension '.h5ad'.")
  SeuratDisk::SaveH5Seurat(SeuratObj, gsub(".h5ad$", ".h5Seurat", h5ad_filename), overwrite = TRUE)
  out_files = c()
  for(assay in assays){
    message(paste0("Writing assay '", assay, "' into .h5ad file."))
    if( ! assay %in% Seurat::Assays( SeuratObj ))
      stop("The assay name given in 'ighc_count_assay_name' is not found in SeuatObj.")
    out_file <- gsub(".h5ad$", paste0("_assay-", assay, ".h5ad"), h5ad_filename)
    SeuratDisk::Convert(gsub(".h5ad$", ".h5Seurat", h5ad_filename), out_file,
                        assay = assay, overwrite = TRUE)
    # see https://github.com/theislab/scvelo/issues/255
    # here read the adata in and perform the fix
    use_condaenv(conda_env, required = TRUE)
    py_run_string("import scanpy as sc")
    py_run_string(paste0("adata = sc.read_h5ad('", out_file, "')"))
    py_run_string("adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})")
    py_run_string(paste0("adata.write_h5ad('", out_file, "')"))
    out_files <- c(out_files, out_file)
  }
  return( out_files )
}

#' parse the substring inside a given cell identifier which corresponds to the nucleotide barcode
#'
#' @description
#' \code{guessBarcodes} parses given cell identifiers to identify the substring which correspond to the nucleotide barcode included in the experiment.
#'
#' @details
#' Numeric / string prefices/suffices were typically added to cell identifiers to avoid wrong mapping across samples; however often these create issues when trying
#' to merge data on the *same* sample but annotated using different workflows.
#' This function attempts to resolve such issues by extracting the nucleotide barcodes actually introduced in the experiment.
#'
#' @param cell_name character, a cell identifier, typicall with prefix and/or suffix (e.g. "ACTGATGCAT-1", "SampleA_ATGAACCTATGG")
#' @param min_barcode_length minimum length of the nucleotide barcode (Default: 6)
#'
#' @return a vector with the input \code{cell_name} decomposed into these three entries:
#' \describe{
#'   \item{prefix}{prefix which exists in the input \code{cell_name} (\code{NA} if doesn't exist in \code{cell_name})}
#'   \item{cell_name}{the actual nucleotide barcode}
#'   \item{suffix}{suffix which exists in the input \code{cell_name} (\code{NA} if doesn't exist in \code{cell_name})}
#' }
#' @importFrom stringr str_detect str_locate str_sub
#' @export guessBarcodes
guessBarcodes <- function(cell_name, min_barcode_length = 6)
{
  # guess the barcode contained in cell_name by looking for continuous
  # nucleotide-like strings in 'cell_name'
  # return a list of c(prefix, cell_name, suffix) (NA if none)
  if( ! stringr::str_detect(cell_name,
                            paste0("[ATCG]{", min_barcode_length, ",}")) ){
    stop("the given cell_name doesn't appear to contain nucleotide barcode strings. Need a custom way to extract cell barcodes.")
  }
  # it would only return the first instance
  barcode_pos <- stringr::str_locate(
    cell_name, paste0("[ATCG]{", min_barcode_length, ",}")
  )
  barcode <- apply(barcode_pos, MARGIN = 1,
                   function(x) stringr::str_sub(cell_name, x[1], x[2]))
  if((barcode_pos[1, 1] - 1) > 1){
    prefix <- stringr::str_sub(cell_name, 1, barcode_pos[1, 1] - 1)
  } else prefix <- NA
  if((barcode_pos[1, 2] + 1) <  nchar(cell_name)){
    suffix <- stringr::str_sub(cell_name, barcode_pos[1, 2] + 1, nchar(cell_name))
  } else suffix <- NA
  return(c(prefix, barcode, suffix))
}

#' Combine velocyto loom files on multiple BAM files into one loom file
#'
#' @description
#' \code{combineLoomFiles} combines .loom files generated using velocyto, on multiple BAM files,
#' into one loom file with the cell barcodes fixed to reflect the cell names in the given Seurat object.
#'
#' @details
#' \code{combineLoomFiles} take a vector of \code{sample_names} (which is assumed to be of the same length and
#' in the same order as \code{loom_files}), parse the prefix and suffix added to the cell barcodes belonging to the
#' given sample, and modify the column names of the matrices in the loom files accordingly. This allows for
#' combining the matrices coming from different samples without ambiguity of cell barcodes.
#' It then checks the matrices for overlap of genes with the given Seurat Object, and remove duplicated genes,
#' and finally write out the merged loom matrices into a new loom file, to be used for RNA velocity analysis.
#'
#' @param loom_files vector of loom files to be merged
#' @param new_loom_filename, character, name of the new loom file to be written out containing the merged data
#' @param SeuratObj corresponding Seurat Object
#' @param sample_names vector of sample names to be looked for in the \code{seurat_sample_column} column of the Seurat metadata.
#' Assumed this is in the same order as the order given in \code{loom_files}.
#' @param seurat_sample_column column name in the Seurat metadata where the different values given in \code{sample_names} can be found
#'
#' @return a output message indicating success of writing out the merged loom matrices into the file given by \code{new_loom_file}.
#'
#' @import hdf5r
#' @export combineLoomFiles
combineLoomFiles <- function(loom_files, new_loom_filename,
                             SeuratObj, sample_names,
                             seurat_sample_column = 'sample_id')
{
  # Here establish a mapping table from the Seurat object and then
  # overwrite the barcodes in the loom objects generated by velocyto CLI
  # to merge them into one combined loom object
  # return new_loom_filename (the merged, barcode-repaired loom file)

  if( length(loom_files) != length(sample_names) )
    stop("'loom_files' and 'sample_names' should be two vectors of the same length, and assumbed to be in the same order.")
  loom_objs <- lapply(loom_files, read_loom_matrices)
  # strip the -[0-9] suffix in the VDJ data frame barcode column
  barcodes_loom <- list()
  message("Assuming the order in sample_names correspond to the order in loom_files. If this is not the case please rerun this function ensuring the order of these vectors match up.")
  for (i in seq_along(sample_names)) {
    barcodes_loom[[i]] <- sapply(colnames(loom_objs[[i]]$spliced), guessBarcodes)
  }
  names( barcodes_loom ) <- sample_names
  barcodes_SeuratObj <- list()
  for (i in seq_along(sample_names)) {
    md <- SeuratObj@meta.data[which(SeuratObj@meta.data[, seurat_sample_column] == sample_names[i]), ]
    barcodes_SeuratObj[[i]] <- sapply(rownames(md), guessBarcodes)
  }
  names( barcodes_SeuratObj ) <- sample_names
  for(i in seq_along(sample_names)){
    loom <- loom_objs[i]
    sample_name <- sample_names[i]
    prefix_Seurat <- unique(barcodes_SeuratObj[[sample_name]][1, ])
    if( length(prefix_Seurat) > 1 ){
      stop(paste0("Sample '", sample_name, "' has more than one prefix in the cell names in SeuratObj. Please fix this."))
    } else if( is.na(prefix_Seurat) ){
      prefix_Seurat <- NULL
    }
    suffix_Seurat <- unique(barcodes_SeuratObj[[sample_name]][3, ])
    if( length(suffix_Seurat) > 1 ){
      stop(paste0("Sample '", sample_name, "' has more than one suffix in the cell names in SeuratObj. Please fix this."))
    } else if( is.na(suffix_Seurat) ){
      suffix_Seurat <- NULL
    }
    barcodes_loom[[sample_name]] <- apply(barcodes_loom[[sample_name]], MARGIN = 2, function(x){
      paste0(prefix_Seurat, x[2], suffix_Seurat)
    })
    colnames(loom_objs[[i]]$spliced) <- barcodes_loom[[sample_name]]
    colnames(loom_objs[[i]]$unspliced) <- barcodes_loom[[sample_name]]
    colnames(loom_objs[[i]]$ambiguous) <- barcodes_loom[[sample_name]]
  }

  # merge loom objects
  new_loom <- list(spliced = as.matrix(Reduce(function(m1, m2) cbind(m1, m2), lapply(loom_objs, function(z) z$spliced))),
                   unspliced = as.matrix(Reduce(function(m1, m2) cbind(m1, m2), lapply(loom_objs, function(z) z$unspliced))),
                   ambiguous = as.matrix(Reduce(function(m1, m2) cbind(m1, m2), lapply(loom_objs, function(z) z$ambiguous))))
  new_loom <- list(
    spliced = new_loom$spliced[ rownames(new_loom$spliced) %in% rownames(SeuratObj@assays$RNA),
                                colnames(new_loom$spliced) %in% colnames(SeuratObj@assays$RNA)],
    unspliced = new_loom$unspliced[ rownames(new_loom$unspliced) %in% rownames(SeuratObj@assays$RNA),
                                    colnames(new_loom$unspliced) %in% colnames(SeuratObj@assays$RNA)],
    ambiguous = new_loom$ambiguous[ rownames(new_loom$ambiguous) %in% rownames(SeuratObj@assays$RNA),
                                    colnames(new_loom$ambiguous) %in% colnames(SeuratObj@assays$RNA)]
  )

  rna_mat <- as.matrix(SeuratObj@assays$RNA@counts)
  # remove all genes which are duplicated in the loom assay due to identical gene names
  # mapped to two distinct Ensembl genes. Removal from both rna_mat and the loom matrices
  duplicated_genes <- rownames(new_loom$spliced)[which(duplicated(rownames(new_loom$spliced)))]
  for(dup_gene in duplicated_genes){
    new_loom$spliced <- new_loom$spliced[which(rownames(new_loom$spliced) != dup_gene), ]
    new_loom$unspliced <- new_loom$unspliced[which(rownames(new_loom$unspliced) != dup_gene), ]
    new_loom$ambiguous <- new_loom$ambiguous[which(rownames(new_loom$ambiguous) != dup_gene), ]
    rna_mat <- rna_mat[which(!grepl(dup_gene, rownames(rna_mat))), ]
  }

  # remove genes which only exist in new_loom / rna_mat
  for(gene in rownames(new_loom$spliced)[which(! rownames(new_loom$spliced) %in% rownames(rna_mat))]){
    new_loom$spliced <- new_loom$spliced[which(rownames(new_loom$spliced) != gene), ]
    new_loom$unspliced <- new_loom$unspliced[which(rownames(new_loom$unspliced) != gene), ]
    new_loom$ambiguous <- new_loom$ambiguous[which(rownames(new_loom$ambiguous) != gene), ]
  }
  for(gene in rownames(rna_mat)[which(! rownames(rna_mat) %in% rownames(new_loom$spliced))]){
    rna_mat <- rna_mat[which(rownames(rna_mat) != gene), ]
  }

  # remove cells from rna_mat which does not exist in new_loom
  rna_mat <- rna_mat[, which(colnames(rna_mat) %in% colnames(new_loom$spliced))]

  # Access loom metadata attributes
  file.orig <- H5File$new(loom_files[1], mode = 'r')

  # write out in hdf5 format as needed by velocyto.R::read.loom.matrices
  file.h5 <- H5File$new(new_loom_filename, mode="w")
  file.h5$create_group("attrs")
  file.h5[["attrs/CreationDate"]] <- file.orig[["attrs"]][["CreationDate"]][]
  file.h5[["attrs/LOOM_SPEC_VERSION"]] <- file.orig[["attrs"]][["LOOM_SPEC_VERSION"]][]
  file.h5[["attrs/velocyto.__version__"]] <- file.orig[["attrs"]][["velocyto.__version__"]][]
  file.h5[["attrs/velocyto.logic"]] <- file.orig[["attrs"]][["velocyto.logic"]][]
  file.h5$create_group("row_attrs")
  file.h5[["row_attrs/Gene"]] <- rownames(new_loom$spliced)
  file.h5$create_group("col_attrs")
  file.h5[["col_attrs/CellID"]] <- colnames(new_loom$spliced)
  file.h5$create_group("layers")
  file.h5[["layers/spliced"]] <- t(new_loom$spliced)
  file.h5[["layers/unspliced"]] <- t(new_loom$unspliced)
  file.h5[["layers/ambiguous"]] <- t(new_loom$ambiguous)
  #file.h5$create_group("matrix")
  file.h5[["matrix"]] <- t(rna_mat)

  ## Close the file at the end
  ## the 'close' method closes only the file-id, but leaves object inside the file open
  ## This may prevent re-opening of the file. 'close_all' closes the file and all objects in it
  file.h5$close_all()
  file.orig$close_all()

  return(paste0("Loom object with velocyto spliced/unspliced counts written to file '",
                new_loom_filename, "'."))
}

#' Combine velocyto loom data with AnnData
#'
#' @description
#' \code{mergeVelocytoWithGEX} merges the velocyto spliced/unspliced gene counts with the AnnData object holding single-cell gene expression data.
#' This is the preprocessing function before calculating RNA velocity using the python scVelo package and workflow.
#'
#' @details
#' \code{mergeVelocytoWithGEX} uses the R reticulate package to run python commands which merge the velocyto counts with the anndata object.
#' It is assumed that the cell barcodes in the velocyto loom and the anndata object matches up. Please consult \code{\link{guessBarcodes}} and \code{\link{combineLoomFiles}} first to ensure these barcodes do match up.
#'
#' @param anndata_file filename pointing to the AnnData file containing gene expression data.
#' @param loom_file, character, name of the new loom file to be written out containing the merged data
#' @param anndata_out_filename output filename of the merged AnnData object to be written containing both the velocyto data and the gene expression data.
#' @param conda_env character, if not \code{NULL} this named conda environment is used to perform the merge.
#' (Default: 'scicsr'). If \code{NULL}, no conda environment will be used, the program assumes the python packages \code{scanpy} and \code{scvelo} are installed in the local python)
#'
#' @return a output message indicating success of writing out the merged AnnData object into the file given by \code{anndata_out_filename}.
#'
#' @import reticulate
#' @export mergeVelocytoWithGEX
mergeVelocytoWithGEX <- function(anndata_file, loom_file, anndata_out_filename,
                                 conda_env = 'scicsr')
{
  use_condaenv(conda_env, required = TRUE)
  py_run_string("import scanpy as sc")
  py_run_string("import scvelo as scv")
  py_run_string(paste0("adata = sc.read_h5ad('", anndata_file, "')"))
  py_run_string("sc.tl.pca(adata)")
  py_run_string("sc.pp.neighbors(adata)")
  message("merging GEX data with velocyto counts ...")
  py_run_string(paste0("adata_velo = scv.read('", loom_file, "')"))
  py_run_string("adata = scv.utils.merge(adata, adata_velo)")
  py_run_string(paste0("adata.write_h5ad('", anndata_out_filename, "')"))
  return(paste0("AnnData file with velocyto counts merged is written to ",
                anndata_out_filename, ""))
}

#' Calculate RNA velocity using the python scvelo workflow
#'
#' @description
#' \code{run_scVelo} calculates RNA velocity using the python scvelo package and the velocity models it implements.
#'
#' @details
#' \code{run_scVelo} uses the R reticulate package to run python commands which run scvelo RNA velocity calculations.
#' It follows the ["RNA Velocity Basics"](https://scvelo.readthedocs.io/VelocityBasics/) tutorial in the scvelo
#' documentation. Unfortunately due to conflicts of the plotting functionalities of R and python this function
#' does **NOT** implement the visualisation of velocity stream onto the dimensionality-reduced projection.
#'
#' @param anndata_file filename pointing to the AnnData file containing gene expression data and merged velocyto spliced/unspliced gene counts.
#' @param anndata_out_filename output filename of the merged AnnData object to be written the fitted RNA velocity estimates calculated using scvelo.
#' @param conda_env character, if not \code{NULL} this named conda environment is used to run scVelo.
#' (Default: 'scicsr'). If \code{NULL}, no conda environment will be used, the program assumes the python packages \code{scanpy} and \code{scvelo} are installed in the local python)
#' @param scvelo_mode the 'mode' parameter in the python scvelo function \code{scv.tl.velocity}. (Default: "dynamical")
#' @param reduction the dimensionality reduction to project RNA velocity estimates onto (Default: "UMAP")
#' @param min_shared_counts include only genes detected in at least this number of cells. (Default: 20)
#' @param n_top_genes include only this many genes with the largest dispersion in the dataset (Default: 2000)
#'
#' @return a output message indicating success of writing out the AnnData object with merged scVelo results into the file given by \code{anndata_out_filename}.
#'
#' @import reticulate
#' @importFrom stringr str_to_lower
#' @export run_scVelo
run_scVelo <- function(anndata_file, anndata_out_filename, conda_env = 'scicsr',
                       scvelo_mode = "dynamical", reduction = "UMAP",
                       min_shared_counts = 20, n_top_genes = 2000)
{
  use_condaenv(conda_env, required = TRUE)  
  arguments <- paste0(
    paste0( system.file( package = "sciCSR" ), "/python/run_scvelo.py" ),
    " --anndata_file ", anndata_file,
    " --anndata_out_file ", anndata_out_filename,
    " --mode ", scvelo_mode,
    " --min_shared_counts ", min_shared_counts,
    " --n_top_genes ", n_top_genes,
    " --reduction ", stringr::str_to_lower(reduction)
  )
  system2(command = py_config()[["python"]], args = arguments)
}

#' Split AnnData object by levels in a specified meta data trait
#'
#' @description
#' \code{splitAnnData} splits a AnnData object by a given metadata column and write out the splitted data subsets into separate .h5ad files.
#'
#' @details
#' \code{splitAnnData} uses the R reticulate package to run python commands which subset the data by the metadata column
#' ( i.e. column in adata.obs) and write out each susbet as separate AnnData objects in .h5ad files.
#'
#' @param anndata_file filename pointing to the AnnData file.
#' @param split.by column name in the metadata (i.e. in \code{SeuratObj@meta.data} considering the Seurat object, or \code{adata.obs} considering the AnnData object)
#' indicating distinct levels to split the AnnData.
#' @param levels vector of values which can be found in the column given by \code{split.by}. Subsets of the Anndata by each value of \code{levels} will be made and written out into .h5ad files.
#' @param conda_env character, if not \code{NULL} this named conda environment is used to perform the split.
#' (Default: 'scicsr'). If \code{NULL}, no conda environment will be used, the program assumes the python packages \code{scanpy} and \code{scvelo} are installed in the local python)
#'
#' @return a vector of output .h5ad filenames with the indicated subset of the data (given in the filenames, see examples.)
#'
#' @import reticulate
#' @export splitAnnData
splitAnnData <- function(anndata_file, split.by, levels, conda_env = 'scicsr')
{
  # split the AnnData object by a given column ('split.by') in Anndata.obs,
  # given the set of values ('levels') to subset for.
  # Write out separate .h5ad file of the subsetted AnnData object
  use_condaenv(conda_env, required = TRUE)
  py_run_string("import scanpy as sc")
  py_run_string("import warnings")
  py_run_string("warnings.filterwarnings('ignore')")
  py_run_string(paste0("adata = sc.read_h5ad('", anndata_file, "')"))
  out_fn <- c()
  for(lvl in levels){
    message(paste0("Subsetting AnnData for ", split.by,  " == '", lvl, "' ..."))
    py_run_string(paste0("adata_subset = adata[adata.obs['", split.by, "'] == '", lvl, "']"))
    subset_filename <- paste0(gsub(".h5ad", "", basename(anndata_file)), "_", lvl, ".h5ad")
    subset_filename <- paste0(dirname(anndata_file), "/", fs::path_sanitize(subset_filename))
    # see https://github.com/theislab/scvelo/issues/255
    py_run_string("adata_subset.__dict__['_raw'].__dict__['_var'] = adata_subset.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})")
    py_run_string(paste0("adata_subset.write_h5ad('", subset_filename, "')"))
    out_fn <- c(out_fn, subset_filename)
  }
  return( out_fn )
}

#' Fit transition model on data using the python cellrank package
#'
#' @description
#' \code{fitTransitionModel} fits transition models on the data using the python cellrank package.
#'
#' @details
#' \code{fitTransitionModel} currently implements either the velocity kernel in cellrank (i.e. uses RNA velocity
#' information to fit transition probabilities) or the pseudotime kernel; a user-indicated column in the metadata
#' will be used as pseudotime reference to fit transition probabiltiies.
#' **NOTE:** In cases where the warning 'Biased KNN graph is disconnected' which will subsequently cause the \code{\link{fitTPT}()} function in the pipeline to falied with error, in our experience it is likely to be caused by subsetting the data prior to computing transitions. Try setting \code{do_pca = FALSE} will preserve the original PCA and avoid this error.
#'
#' @param anndata_file filename pointing to the AnnData file.
#' @param conda_env character, if not \code{NULL} this named conda environment is used to run CellRank.
#' (Default: \code{NULL}, i.e. no conda environment will be used, the program assumes the python packages \code{scanpy}, \code{scvelo} and \code{cellrank} are installed in the local python)
#' @param mode character, either 'pseudotime' (uses the cellrank 'PseudotimeKernel') or 'velocity' (cellrank 'VelocityKernel'). (Default: 'pseudotime')
#' @param pseudotime_key character, column name which indicates the ranking to be used as pseudotime ordering of the cells. Not considered if mode is 'velocity'. (Default: 'csr_pot')
#' @param do_pca Should principal component analysis (PCA) be re-computed on the data? (Default: TRUE)
#' @param do_neighbors Should k-nearest neighbour (kNN) graph be re-computed on the data? (Default: TRUE)
#'
#' @return a list with three entries:
#' \describe{
#'   \item{cellrank_obj}{Python \code{cellrank.tl.estimators.CFLARE} object containing details of the fitted transition model}
#'   \item{transition_matrix}{The matrix holding cell-to-cell transition probability (i.e. \code{cellrank_obj.transition_matrix}), but converted into a dense matrix. This is the single-cell transition matrix for downstream uses.}
#'   \item{CellID}{a vector of cell identifiers in the order of each row/column of \code{transition_matrix}.}
#' }
#'
#' @import reticulate
#' @export fitTransitionModel
fitTransitionModel <- function(anndata_file, conda_env = 'scicsr', mode = 'pseudotime',
                               pseudotime_key = 'csr_pot', do_pca = TRUE, do_neighbors = TRUE)
{
  if( ! mode %in% c("pseudotime", "velocity"))
    stop("Currently 'mode' must be either 'pseudotime' (cellrank PseudotimeKernel) or 'velocity' (cellrank VelocityKernel).")
  use_condaenv(conda_env, required = TRUE)
  # update executable path in sys module if OS is windows
  # see https://github.com/rstudio/reticulate/issues/517
  if( .Platform$OS.type == "windows" ){
    sys <- import("sys")
    exe <- file.path(sys$exec_prefix, "pythonw.exe")
    sys$executable <- exe
    sys$`_base_executable` <- exe
    
    # update executable path in multiprocessing module
    multiprocessing <- import("multiprocessing")
    multiprocessing$set_executable(exe)
  }
  
  py_run_string("import scanpy as sc")
  py_run_string(paste0("adata = sc.read_h5ad('", anndata_file, "')"))
  if( do_pca ){
    py_run_string("sc.tl.pca(adata)")
  }
  if( do_neighbors ){
    py_run_string("sc.pp.neighbors(adata)")
  }
  fit_cellrank_pseudotime_kernel <- NULL
  fit_cellrank_velocity_kernel <- NULL
  source_python(paste0( system.file( package = "sciCSR" ), "/python/cellrank_functions.py" ) )
  if( mode == 'pseudotime' ){
    g <- fit_cellrank_pseudotime_kernel(py$adata, pseudotime_key = pseudotime_key)
  } else if (mode == "velocity" ){
    g <- fit_cellrank_velocity_kernel(py$adata)
  }
  return( g )
}

#' Fit Transition Path Theory (TPT) on the cellrank transition models
#'
#' @description
#' \code{fitTPT} fits Transition Path Theory on the transition model defined using \code{\link{fitTransitionModel}()}, at
#' a 'coarse-grained' level where transitions are considered between *groups* of cells with grouping indicated by the user.
#'
#' @details
#' \code{fitTPT} interfaces with (and reimplements some routines to improve efficincy) the Python \code{deeptime} package to fit transition path theory (TPT) onto the markov state model defined by running
#' the \code{\link{fitTransitionModel}} function that uses cellrank under the hood. With the parameter \code{group.cells.by}, the user specifies
#' a scheme to group individual row/columns of the transition matrix (for example, by cell type or by isotype). The function
#' then fits TPT on to this grouped/'coarse-grained' transition matrix, upon user indicating a likely 'source' and 'target' state.
#' The output are estimated information flows ('flux') between different states in order to flow from the source to the target,
#' and the probabilities of sampling each state at equilibrium ('stationary distribution').
#' A random 'null background' model was fitted by randomly reshuffling columns of the transition matrix by \code{random_n} (default: 100) times.
#' These random fluxes help determine the significance of an observed flux, by calculating one-sided empirical probabilities of the observed flux larger than the that observed in the randomised models.
#'
#' @param anndata_file filename pointing to the AnnData file.
#' @param CellrankObj the \code{cellrank_obj} entry of the output list of \code{\link{fitTransitionModel}()}.
#' @param group.cells.by character, column in the metadata to group cells
#' @param source_state character, a value in the \code{group.cells.by} column which is taken as the source state for fitting transition path theory. All cells belonging to this group are considered as the source.
#' @param target_state character, a value in the \code{group.cells.by} column which is taken as the target state for fitting transition path theory. All cells belonging to this group are considered as the target.
#' @param conda_env character, if not \code{NULL} this named conda environment is used to perform TPT analysis.
#' (Default: \code{NULL}, i.e. no conda environment will be used, the program assumes the python packages \code{scanpy}, \code{scvelo} and \code{cellrank} are installed in the local python)
#' @param random_n number of times to reshuffle transition matrix columns to derive randomised models (default: 100).
#' @param do_pca Should principal component analysis (PCA) be re-computed on the data? (Default: TRUE)
#' @param do_neighbors Should k-nearest neighbour (kNN) graph be re-computed on the data? (Default: TRUE)
#'
#' @return a list with these entries:
#' \describe{
#'   \item{gross_flux}{a n-by-n matrix (where n is the total number of states), of total fluxes estimated between from a state (row) to another state (column). }
#'   \item{pathways}{a data.frame indicating the possible paths to take from \code{source_state} to \code{target_state}, and the likelihood (max: 100) to travel through each stated path.}
#'   \item{significance}{a n-by-n matrix (where n is the total number of states), where the observed gross flux is greater than the flux estimated in the randomised models. }
#'   \item{total_gross_flux}{element-wise sum of the gross_flux matrix. }
#'   \item{total_gross_flux_reshuffled}{element-wise sum of the gross_flux matrix, calculated over each randomised (randomly reshuffled transition matrix coluns) models.}
#'   \item{gross_flux_randomised}{gross_flux matrix but from the randomised (randomly reshuffled transition matrix coluns) TPT models.}
#'   \item{mfpt}{Mean First Passage Time required to travel from \code{source_state} to \code{target_state} as estimated by Transition Path Theory. }
#'   \item{mfpt_reshuffled}{Mean First Passage Time required to travel from \code{source_state} to \code{target_state} as estimated by Transition Path Theory, calculated over each randomised (randomly reshuffled transition matrix coluns) models. }
#'   \item{stationary_distribution}{Equilibrium probability of each state as estimated by Transition Path Theory. }
#'   \item{stationary_distribution_reshuffled}{Equilibrium probability of each state as estimated by Transition Path Theory, calculated over each randomised (randomly reshuffled transition matrix coluns) models. }
#' }
#'
#' @import reticulate
#' @export fitTPT
fitTPT <- function(anndata_file, CellrankObj,
                   group.cells.by, source_state, target_state,
                   conda_env = 'scicsr',random_n = 100,
                   do_pca = TRUE, do_neighbors = TRUE)
{
  fit_coarse_grain_tpt <- NULL
  use_condaenv(conda_env, required = TRUE)
  py_run_string("import scanpy as sc")
  source_python(paste0( system.file( package = "sciCSR" ), "/python/TPT_functions.py" ) )
  py_run_string(paste0("adata = sc.read_h5ad('", anndata_file, "')"))
  py_run_string(paste0("cluster_ident = adata.obs.groupby('",
                       group.cells.by, "').indices"))
  py_run_string("cluster_names = [i for i in cluster_ident.keys()]")
  cluster_names <- sapply(py$cluster_names, as.character)
  tpt <- fit_coarse_grain_tpt(CellrankObj$cellrank_obj$transition_matrix,
                              py$cluster_ident,
                              source_state, target_state, random_n = random_n)
  if( ! "data.frame" %in% class(tpt[["pathways"]]) ){
    tpt[["pathways"]] <- reticulate::py_to_r(tpt[["pathways"]])
  }
  tpt[["stationary_distribution_bootstrapping"]] <- lapply(tpt[["stationary_bootstrapping"]], unlist)
  tpt[["total_gross_flux"]] <- tpt[["total_gross_flux"]]
  tpt[["total_gross_flux_reshuffled"]] <- unlist(tpt[["total_gross_flux_randomised"]])
  tpt[["gross_flux_randomised"]] <- lapply(1:random_n, function(i) tpt[["gross_flux_randomised"]][i, , ])
  tpt[["mfpt"]] <- tpt$coarse_grain_tpt$mfpt
  tpt[["mfpt_reshuffled"]] <- sapply(tpt$randomised_tpt, function(x) x$tpt$mfpt)

  # name the clusters and reorder them in the order of the
  # original factor indicating the grouping
  cluster_order <- c(source_state, cluster_names[! cluster_names %in% c(source_state, target_state)],
                     target_state)
  tpt[["stationary_distribution"]] <- c(tpt$coarse_grain_tpt$stationary_distribution)
  names(tpt[["stationary_distribution"]]) <- cluster_order
  tpt[["stationary_distribution"]] <- tpt[["stationary_distribution"]][cluster_names]
  tpt[["stationary_distribution_reshuffled"]] <- lapply(tpt$randomised_tpt, function(x){
    o <- c(x$tpt$stationary_distribution)
    names(o) <- cluster_order
    o[cluster_names]
  })
  dimnames(tpt[["gross_flux"]]) <- list(cluster_order, cluster_order)
  tpt[["gross_flux"]] <- tpt[["gross_flux"]][cluster_names, cluster_names]
  dimnames(tpt[["significance"]]) <- list(cluster_order, cluster_order)
  tpt[["significance"]] <- tpt[["significance"]][cluster_names, cluster_names]
  for(i in 1:random_n){
    dimnames(tpt[["gross_flux_randomised"]][[i]]) <- list(cluster_order, cluster_order)
    tpt[["gross_flux_randomised"]][[i]] <-tpt[["gross_flux_randomised"]][[i]][cluster_names, cluster_names]
  }
  return(tpt[c("gross_flux", "pathways", "significance",
               "total_gross_flux", "total_gross_flux_reshuffled",
               "gross_flux_randomised", "mfpt", "stationary_distribution",
               "mfpt_reshuffled", "stationary_distribution_reshuffled",
               "stationary_distribution_bootstrapping")])
}

#' Visualise flux matrix describing class-switch recombination (CSR) transitions in data
#'
#' @description
#' \code{plotFluxMatrix} parse data from \code{\link{fitTPT}()} and visualises the estimated isotype-switching dynamics in the form of a flux matrix detailing amount of switches from and to each isotype.
#'
#' @details
#' \code{plotFluxMatrix} parses data from the \code{\link{fitTPT}()} function and visualises the results from Transition Path Theory (TPT) in the form of a bubble plot representing the amount of flux from and to each isotype.
#' * The bubble size is scaled by significance, i.e. the likelihood that the observed flux is greater than randomised flux estimates obtained by reshuffling columns of the transition matrix.
#' * The bubble colour is scaled by magnitude, i.e. the amount of flux estimated to flow from one isotype to another isotype.
#' * Improbable CSR combinations (i.e. switching back to an isotype preceding the current isotype) are removed by default. This can be turned off by indicating \code{mask_improbable_csr = FALSE}.
#'
#' @param TPTObj List of TPT results, output from the \code{\link{fitTPT}} function.
#' @param SeuratObj Seurat object
#' @param ighc_count_assay_name name of assay in SeuratObj which holds the IgH productive/sterile transcript counts. (Default: "IGHC")
#' @param mask_improbable_csr Should isotype combinations which represents improbable Class-switch recombination events (i.e. switching back to an isotype 5' to the current one) be removed from visualisation? (Default: TRUE)
#' @param return_plot Should the CSR transition plot be returned? If FALSE, a named list of \code{stationary_distribution} and \code{flux} will be returned which contains the data frames to be visualised in this plot. (Default: TRUE)
#'
#' @return A ggplot2 object showing the flux matrix and the associated significance level in the form a bubble plot. See 'Description' for details.
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom reshape2 melt
#' @importFrom stringr str_to_lower str_replace str_to_sentence
#' @importFrom stats p.adjust
#' @import ggplot2
#' @export plotFluxMatrix
plotFluxMatrix <- function(TPTObj, SeuratObj,
                           ighc_count_assay_name = "IGHC",
                           mask_improbable_csr = TRUE,
                           return_plot = TRUE)
{
  flux_matrix <- TPTObj$gross_flux
  significance_matrix <- TPTObj$significance

  # check whether the flux matrix contains all isotypes; if not, add them back
  ighc_counts <- SeuratObj@assays[[ighc_count_assay_name]]@counts
  c_genes <- rownames(ighc_counts)[grepl("-C$", rownames(ighc_counts))]
  c_genes <- stringr::str_replace(stringr::str_replace(c_genes, "-C$", ""),
                                  "^IGH|^Igh", "")
  c_genes <- stringr::str_to_sentence(c_genes)
  # check whether row/column names of flux_matrix conform to this format
  # if not , force them to by removing IGH or Igh
  if(sum(sapply(c_genes, function(c_gene) c_gene %in% rownames(flux_matrix))) == 0){
    rownames(flux_matrix) <- stringr::str_replace(rownames(flux_matrix), "^IGH|^Igh", "")
    rownames(flux_matrix) <- stringr::str_to_sentence(stringr::str_to_lower(rownames(flux_matrix)))
    colnames(flux_matrix) <- stringr::str_replace(colnames(flux_matrix), "^IGH|^Igh", "")
    colnames(flux_matrix) <- stringr::str_to_sentence(stringr::str_to_lower(colnames(flux_matrix)))
  }
  for(c_gene in c_genes){
    if( ! c_gene %in% rownames(flux_matrix) ){
      flux_matrix <- rbind(flux_matrix, rep(0, ncol(flux_matrix)))
      rownames(flux_matrix)[nrow(flux_matrix)] <- c_gene
      flux_matrix <- cbind(flux_matrix, rep(0, nrow(flux_matrix)))
      colnames(flux_matrix)[ncol(flux_matrix)] <- c_gene
      significance_matrix <- rbind(significance_matrix, rep(0, ncol(significance_matrix)))
      rownames(significance_matrix)[nrow(significance_matrix)] <- c_gene
      significance_matrix <- cbind(significance_matrix, rep(0, nrow(significance_matrix)))
      colnames(significance_matrix)[ncol(significance_matrix)] <- c_gene
    }
  }
  flux_matrix <- flux_matrix[c_genes, c_genes]
  significance_matrix <- significance_matrix[c_genes, c_genes]
  # filtering
  #  flux_matrix[significance_matrix < significance_threshold] <- 0
  if(mask_improbable_csr){
    flux_matrix[lower.tri(flux_matrix)] <- 0
  }
  #  if( is.numeric(mask_threshold) ) {
  #    flux_matrix[flux_matrix <= mask_threshold] <- 0
  #  }
  graph <- reshape2::melt(flux_matrix, value.name = "flux", varnames = c("from", "to"))
  graph <- graph[which(graph$flux > 0), ]
  graph <- merge(graph,
                 reshape2::melt(significance_matrix, varnames = c("from", "to"),
                                value.name = "pval"))
  graph$signif <- -log(p.adjust(graph$pval))
  p <- ggplot(graph, aes_string(x = "from", y = "to", color = "flux", size = "signif")) +
    geom_point() + cowplot::theme_cowplot() + scale_colour_gradient2(name = "% flux") +
    scale_size_continuous(name = "p-value", breaks = c(-log(1), -log(0.5), -log(0.1), -log(0.05)),
                          labels = c("1", "0.5", "0.1", "0.05")) +
    scale_x_discrete(drop = FALSE, position = "top") + scale_y_discrete(drop = FALSE)
  if( return_plot ) return(p)
  else return(graph)
}

#' Visualise stationary distribution of isotypes in the data
#'
#' @description
#' \code{plotStationaryDistribution} parse data from \code{\link{fitTPT}()} to plot a bar-plot of isotypes in the given data. The stationary distribution represents the equilibrium distribution of cells harbouring each isotype, taking the underlying isotype-switching dynamics into consideration.
#'
#' @details
#' \code{plotStationaryDistribution} parses data from the \code{\link{fitTPT}()} function for visualising the stationary distribution of isotypes. The resulting plot is a bar plot of stationary distribution, with error-bars obtained from bootstrapping (i.e. sampling with replacement) cells harbouring each isotype. The error-bar shown are 95% confidence intervals resulting from the bootstrap sampling.
#'
#' @param TPTObj List of TPT results, output from the \code{\link{fitTPT}()} function.
#' @param SeuratObj Seurat object
#' @param ighc_count_assay_name name of assay in SeuratObj which holds the IgH productive/sterile transcript counts. (Default: "IGHC")
#' @param return_plot Should the CSR transition plot be returned? If FALSE, a named list of \code{stationary_distribution} and \code{flux} will be returned which contains the data frames to be visualised in this plot. (Default: TRUE)
#'
#' @return A ggplot2 object of the bar-plot of stationary distribution with 95% bootstrapped confidence intervals.
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom reshape2 melt
#' @importFrom stringr str_to_lower str_replace str_to_sentence
#' @importFrom stats quantile
#' @import ggplot2
#' @export plotStationaryDistribution
plotStationaryDistribution <- function(TPTObj, SeuratObj,
                                       ighc_count_assay_name = "IGHC",
                                       return_plot = TRUE)
{
  stationary_distribution <- TPTObj$stationary_distribution
  bs_stationary <- TPTObj$stationary_distribution_bootstrapping

  # check whether the flux matrix contains all isotypes; if not, add them back
  ighc_counts <- SeuratObj@assays[[ighc_count_assay_name]]@counts
  c_genes <- rownames(ighc_counts)[grepl("-C$", rownames(ighc_counts))]
  c_genes <- stringr::str_replace(stringr::str_replace(c_genes, "-C$", ""),
                                  "^IGH|^Igh", "")
  c_genes <- stringr::str_to_sentence(c_genes)
  for(c_gene in c_genes){
    if( ! c_gene %in% names(stationary_distribution) ){
      stationary_distribution <- c(stationary_distribution, 0)
      names(stationary_distribution)[length(stationary_distribution)] <- c_gene
    }
  }
  # initialise node annotations
  stationary_distribution <- stationary_distribution[c_genes]
  isotypes <- data.frame(isotype = c_genes, x = 1:length(c_genes), y = 1)
  # calculate bootstrapped 95% confidence intervals
  bs_stationary <- data.frame(t(sapply(bs_stationary, quantile, probs = c(0.025, 0.975))))
  bs_stationary$isotype <- rownames(bs_stationary)
  colnames(bs_stationary) <- c("lowq", "highq", "isotype")
  isotypes <- merge(
    isotypes,
    data.frame("stationary" = stationary_distribution),
    by.x = "isotype", by.y = "row.names", all.x = TRUE, all.y = FALSE, sort = FALSE
  )
  isotypes <- merge(isotypes, bs_stationary,
                    by = "isotype", all.x = TRUE, all.y = FALSE, sort = FALSE)
  isotypes[is.na(isotypes[, "stationary"]), "stationary"] <- 0
  isotypes[is.na(isotypes[, "lowq"]), "lowq"] <- 0
  isotypes[is.na(isotypes[, "highq"]), "highq"] <- 0
  isotypes[, "isotype"] <- factor(isotypes[, "isotype"], levels = c_genes)

  # now do plotting
  p <- ggplot(isotypes) +
    geom_bar(aes_string(x = "isotype", y = "stationary"),
             stat = "identity", position = position_dodge2()) +
    geom_errorbar(aes_string(x = "isotype", ymin = "lowq", ymax = "highq"),
                  width = 0, position = position_dodge2())
  p <- p + scale_y_reverse(name = "stationary\ndistribution") +
    scale_x_discrete(drop=FALSE, position = "top") +
    cowplot::theme_cowplot() #+ scale_x_discrete(breaks = rev(c_genes)) +
  if( return_plot ) return(p)
  else return(isotypes)
}

#' Computing distances between multiple transition matrices
#'
#' @description
#' \code{compareTransitionMatrices} computes distances between a list of transition matrices.
#'
#' @details
#' \code{compareTransitionMatrices} takes the list of transition matrices given in \code{matrix_list}, and sample realisations from
#' the Markov chain defined using the transition matrix. It then compares the similarity of these sampled trajectories
#' using either the Kullback-Leibler divergence (\code{distance_metric == "KL"}) or the Jensen-Shannon divergence (\code{distance_metric == "JSD"}).
#' Both can be interpreted the same way (the larger this number, the more different two transition matrices and their
#' realisations are), although the Jensen-Shannon divergence is scaled between 0 and 1.
#'
#' @param matrix_list list of transition matrices to be compared.
#' @param SeuratObj Seurat object. Considered only if \code{csr_transitions} is \code{NULL}.
#' @param cells vector of cell identifiers corresponding to the row/column order in the transition matrices. (Assumed that all supplied transition matrices have exactly the same row/column ordering.)
#' @param group.by column in the Seurat object metadata on which cells are grouped. (Default: 'seurat_clusters')
#' @param n_realisation number of trajectories ('realisations') to be sampled from the Markov model defined using each transition matrix. (Default: 1000)
#' @param n_step number of time-steps in each trajectory/realisation to be sampled. (Default: 1000)
#' @param distance_metric the distance metric to be calculated. Either "KL" (for Kullback-Leibler divergence) or "JSD" (Jensen-Shannon divergence). (Default: "KL")
#'
#' @return A list of two entries:
#' \describe{
#'   \item{distance}{distance (\code{distance_metric}) between the trajectories sampled from each possible pair of transition matrix.}
#'   \item{sampled_transitions}{list of trajectories sampled from the transition matrices.}
#' }
#' @import markovchain
#' @importFrom methods new
#' @importFrom philentropy KL JSD
#' @export compareTransitionMatrices
compareTransitionMatrices <- function(matrix_list, SeuratObj,
                                      cells, group.by = "seurat_clusters",
                                      n_realisation = 1000, n_step = 1000,
                                      distance_metric = "KL")
{
  if( ! distance_metric %in% c("KL", "JSD") )
    stop("'distance_metric' must be one of 'KL' (Kullback-Leibler divergence, default) or 'JSD' (Jensen-Shannon divergence).")
  mc_list <- list()
  for(i in 1:length(matrix_list)){
    mc_list <- c(mc_list,
                 new("markovchain", states = cells,
                     byrow = TRUE, transitionMatrix = matrix_list[[i]],
                     name = paste0("matrix", i)))
  }

  if( is.factor(SeuratObj[[group.by]][, 1]) ){
    possible_states = levels( SeuratObj[[group.by]][, 1] )
  } else{
    possible_states = sort( unique( SeuratObj[[group.by]][, 1] ) )
  }
  sampled_transitions <- lapply(1:length(mc_list), function(y){
    message(paste0("Sampling realisations from Markov chain ", y, " ..."))
    o <- Reduce("+", lapply(1:n_realisation, function(x){
      set.seed(1000 + x)
      createSequenceMatrix(
        unname(SeuratObj[[group.by]][markovchainSequence(n_step, mc_list[[y]]), 1]),
        possibleStates = possible_states
      )
    }))
    # normalise so that sum of all entries in each sampled_transition matrix = 1
    o <- apply(o, c(1, 2), function(x) x / sum(o))
    return( o )
  })
  names(sampled_transitions) <- paste0("mc_", 1:length(sampled_transitions))
  message("Comparing the list of transition models ...")
  if( distance_metric == "KL" ){
    dist_mat <- philentropy::KL(do.call("rbind", lapply(sampled_transitions, c)))
  } else if( distance_metric == "JSD" ){
    dist_mat <- philentropy::JSD(do.call("rbind", lapply(sampled_transitions, c)))
  }
  if( is.array( dist_mat ) ){
    dimnames(dist_mat) <- list(names(sampled_transitions), names(sampled_transitions))
  }
  return(list("distance" = dist_mat,
              "sampled_transitions" = sampled_transitions))
}

#' Plot arrows showing transitions in the style of RNA velocity plots
#'
#' @description
#' \code{plot_arrows} plots arrows on top of UMAP projection to show transitions in the style of conventional RNA velocity analysis. You can indicate whether to use RNA velocity (\code{based_on = 'velocity'}), CSR (\code{based_on = 'CSR'}) or SHM (\code{based_on = 'SHM'}) information to project the arrows.
#'
#' @details
#' \code{plot_arrows} uses the plotting functionalities in scVelo in python to generate a plot of arrows indicating transitions, given the type of biological information (velocity/CSR/SHM). You can choose between projecting transitions as arrows laid out on a grid (\code{style = 'grid'}) or as streams of arrows (\code{style = 'stream'}). The plot is saved as a SVG/PDF/PNG (depending on file extension given in \code{img_path}, see below), and re-rendered in the 'plot' panel in R.
#'
#' @param anndata_file input anndata_file. If \code{based_on} is 'velocity', this file needs to be output from \code{\link{run_scVelo}}. If \code{based_on} is 'CSR' or 'SHM', the columns 'csr_pot' or 'shm' should be in the .obs slot of the AnnData object.
#' @param img_path Optional, path to write out the arrow plot. If supplied, you can specify file format (PNG/SVG/PDF) by including the file extension. If PNG, an image of 600 dots per inch will be rendered. Default is \code{NULL}, i.e. it will write to a temporary file as a PNG.
#' @param based_on one of 'velocity', 'csr', 'shm'. The type of information to be used to project arrows. Each has requirements on the input \code{anndata_file} (see argument \code{anndata_file} of this function).
#' @param style one of 'grid' (lay out arrows on a grid of fixed width/height on the UMAP plot) or 'stream' (draw arrows as streams), in the style of the scvelo 'pl.velocity_embedding_grid' or 'pl.velocity_embedding_stream' respectively.
#' @param title plot title (Default: \code{NULL}, the title of the plot will be identical to 'based_on')
#' @param colour.by column in gene expression metadata to group and colour the cells by (Default: 'seurat_clusters')
#' @param cols Optional, a vector of characters containing the HEX code of colours to be used. Has to be the same length as the number of levels found in the \code{colour.by} variable in the gene expression metadata.
#' @param components component of dimensionality reduction to show in the plot. e.g. put \'1,3\' if desired UMAP plot displays the UMAP_1 and UMAP_3 axes.(Default: '1,2')
#' @param conda_env character, if not \code{NULL} this named conda environment is used to generate the plot in scVelo.
#' (Default: 'scicsr'). If \code{NULL}, no conda environment will be used, the program assumes the python packages \code{scanpy} and \code{scvelo} are installed in the local python)
#'
#' @return plot of UMAP dimensionality reduction with arrows projected on top depicting inferred transitions. The same plot is saved in the path given by \code{img_path}. If PNG, it is rendered at 600 dots-per-inch (dpi).
#'
#' @importFrom png readPNG
#' @importFrom reticulate py_config use_condaenv
#' @importFrom grid grid.newpage grid.raster
#' @export plot_arrows
plot_arrows <- function(anndata_file, img_path = NULL, based_on = 'velocity',
                        style = "grid", title = NULL,
                        colour.by = 'seurat_clusters',
                        cols = NULL, components = '1,2',
                        conda_env = 'scicsr')
{
  use_condaenv(conda_env, required = TRUE)
  if( is.null( img_path) ){
    img_path <- paste0(tempfile(), ".png")
  }
  arguments <- paste0(
    paste0( system.file( package = "sciCSR" ), "/python/plot_velocity.py" ),
    " --anndata ", anndata_file,
    " --output ", img_path,
    " --type ", based_on,
    " --color ", colour.by,
    " --components ", components,
    " --style ", style
  )
  if( !is.null( cols ) ){
    cols <- paste(cols, collapse = ",")
    arguments <- paste0(arguments, " --palette \"", cols, "\"")
  }
  if( !is.null( title ) ){
    arguments <- paste0(arguments, " --title \"", title, "\"")
  }
  system2(command = py_config()[["python"]], args = arguments)
  img <- readPNG(img_path)
  grid::grid.newpage()
  grid.raster(img)
}
