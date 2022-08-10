#' get somatic hypermutation level
#'
#' @description
#' `getSHM` calculates somatic hypermutation (i.e. 1 - (percentage identity to germline VH gene)).
#'
#' @details
#' `getSHM` considers the `v_identity` (i.e. % identity to germline V gene) for the VH sequence
#' mapped to each cell, and calculate 1 - `v_identity` as the somatic hypermutation (SHM)
#' level for the cell. It finds the v_identity information from the meta.data slot of the Seurat object.
#' For cells without a mapped V sequence, it will impute SHM as 0.
#'
#' @param SeuratObj Seurat Object
#' @param v_identity_anno_name column name in the Seurat meta.data slot which holds the v_identity information to be considered.
#' @param shm_column_to_add name of column to be added to the Seurat meta.data slot
#' which holds the SHM frequency calculated in this function
#'
#' @return The same Seurat object as given by SeuratObj, except that a new column with
#' name given by `shm_column_to_add` is appended to the SeuratObj meta.data slot to reflects
#' the calculated SHM frequency.
#'
#' @importFrom Seurat AddMetaData
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
#' `getCSRpotential` scores each cell by their status in temrs of class switch recombination (CSR), by considering the mapped sterile and productive IgH transcripts..
#'
#' @details
#' `getCSRpotential` calculates a "CSR potential" score which ranks the cells in the given Seurat Obj by their status in the CSR process.
#' This is given by the Euclidean norm of (furthest_jc, total_ic) (i.e.\eqn{ \sqrt{ \text{furthest\_jc}^2 + \text{total\_ic}^2} } ), where
#' \itemize{
#'   \item{furthest_jc}{: the productive isotype for each cell furthest along the IGH loci (0 = IgM, 1 = IgG3 ... ), and}
#'   \item{total_ic}{: amount of sterile IgH molecules for each cell.}
#' }
#' For `total_ic`, the default is to use the scale.data slot which already normalises the IGHC counts by library size. If this doesn't exist the function will calculate this while regressing out the library size.
#' For `furthest_jc`, users can either use a specified column in the Seurat object meta.data which indicates the isotype of the cell, or, if not provided, used the productive reads counted using the productive/sterile quantification workflow implemented in this package.
#'
#' @param SeuratObj Seurat Object
#' @param ighc_count_assay_name name of assay in `SeuratObj` which holds the IgH productive/sterile transcript count data. (Default: "IGHC")
#' @param ighc_slot the slot in `slot(SeuratObj, "assays")[[ighc_count_assay_name]]` to be used to access productive/sterile transcript counts (Default: "scale_data")
#' @param vars.to.regress list of variables to be regressed out in calculating the scale.data slot, if `ighc_slot` is given as `scale.data` but it has not been populated. (Default: "nCount_RNA", i.e. per-cell library size)
#' @param c_gene_anno_name If not NULL, this column from the Seurat Object meta.data will be used to indicate `furthest_jc` in calculaing the CSR potential score, in lieu of the productive transcript counts in the IGHC assay (Default: NULL)
#' @param isotype_column_to_add name of column to be added to the SeuratObj meta.data to indicate the isotype of the cell. Used for subsequent grouping of cells in calculating transitions.
#'
#' @return Seurat object with these following columns added to the meta.data slot:
#' \itemize{
#'   \item{furthest_jc}{an integer indicating the productive isotype for each cell furthest along the IGH loci (0 = IgM, 1 = IgG3 ... )}
#'   \item{total_ic}{amount of sterile IgH molecules for each cell, calculated from the given `ighc_slot` of the IGHC assay.}
#'   \item{csr_pot}{Euclidean norm of furhest_jc and total_ic (i.e. \eqn{ \sqrt{ \text{furthest\_jc}^2 + \text{total\_ic}^2} }), and normalised into range [0, 1] across the given dataset with 0 indicating that the cell at the earliest point in terms of CSR across the dataset.}
#'   \item{`isotype_column_to_add`}{isotype labelled as M, G3, etc. (added only when c_gene_anno_name is FALSE and the ighc_count_assay_name Assay is used to calculate CSR potential.}
#' }
#'
#' @importFrom Seurat Assays ScaleData AddMetaData
#' @importFrom Matrix colSums
#' @importFrom stringr str_replace str_to_sentence str_to_upper
#'
getCSRpotential <- function(SeuratObj, ighc_count_assay_name = "IGHC",
                            ighc_slot = "scale.data",
                            vars.to.regress = c("nCount_RNA"),
                            c_gene_anno_name = NULL,
                            isotype_column_to_add = "isotype")
{
  if(! ighc_count_assay_name %in% Seurat::Assays(SeuratObj))
    stop(paste0("The assay '", ighc_count_assay_name, "' cannot be found in SeuatObj.") )
  if( ! ighc_slot %in% slotNames(SeuratObj@assays[[ighc_count_assay_name]]) )
    stop(paste0("The named slot '", ighc_slot, "' cannot be found in the '",
                ighc_count_assay_name, "' assay in SeuratObj."))
  if( ighc_slot == "scale.data" &&
      all(dim(SeuratObj@assays[[ighc_count_assay_name]]@scale.data) == 0)){
    # populate the scale.data matrix by running Seurat::ScaleData
    SeuratObj <- Seurat::ScaleData(SeuratObj, vars.to.regress = vars.to.regress,
                                   assay = ighc_count_assay_name)
  }
  ighc_counts <- slot(SeuratObj@assays[[ighc_count_assay_name]], ighc_slot)
  ic_genes <- rownames(ighc_counts)[grepl("IC$", rownames(ighc_counts))]
  total_ic <- Matrix::colSums(ighc_counts[ic_genes, ])
  c_genes <- rownames(ighc_counts)[grepl("-C$", rownames(ighc_counts))]
  c_genes <- stringr::str_replace(stringr::str_replace(c_genes, "-C$", ""),
                                  "^IGH|^Igh", "")
  c_genes <- stringr::str_to_sentence(c_genes)

  if( is.null( c_gene_anno_name ) ){
    if( is.null( isotype_column_to_add ) )
      stop("'isotype_column_to_add' needs to be given if 'c_gene_anno_name' is NULL.")
    jc_genes <- rownames(ighc_counts)[grepl("JC$", rownames(ighc_counts))]
    jc_counts <- ighc_counts[jc_genes, ]
    furthest_jc <- apply(jc_counts, MARGIN = 2, function(x) {
      pos <- which(x > 0) # positive entries
      if(length(pos) == 0) return(0)
      else return(max(pos) - 1)
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
    jc_anno <- stringr::str_to_upper(jc_anno)
    c_genes <- stringr::str_replace(stringr::str_replace(ic_genes, "-IC$", ""),
                                    "^IGH|^Igh", "")
    jc_anno <- factor(jc_anno, levels = stringr::str_to_upper(c_genes))
    furthest_jc <- as.numeric(jc_anno) - 1
    furthest_jc[is.na(furthest_jc)] <- 0
    c_gene_anno <- factor(furthest_jc, levels = (1:length(c_genes)) - 1,
                          labels = stringr::str_to_sentence(c_genes))
    names(c_gene_anno) <- colnames(ighc_counts)
    SeuratObj <- Seurat::AddMetaData(SeuratObj, c_gene_anno,
                                     col.name = isotype_column_to_add)
  }
  csr_pot <- sqrt( furthest_jc^2 + total_ic^2)
  csr_pot <- (csr_pot - min(csr_pot, na.rm = TRUE)) / abs(diff(range(csr_pot, na.rm = TRUE)))
  furthest_jc[is.na(furthest_jc)] <- 0
  total_ic[is.na(total_ic)] <- 0
  csr_pot[is.na(csr_pot)] <- 0
  o <- data.frame(furthest_jc = furthest_jc, total_ic = total_ic,
                  csr_pot = csr_pot)
  rownames(o) <- colnames(ighc_counts)
  SeuratObj <- Seurat::AddMetaData(SeuratObj, o)
  return(SeuratObj)
}

#' wrapper function to convert Seurat Object to a AnnData .h5ad file
#'
#' @description
#' `convertSeuratToH5ad` is a wrapper function to convert a given Seurat Object into an AnnData object (for use in python with e.g. scanpy) and write out into a .h5ad file.
#'
#' @details
#' `convertSeuratToH5ad` simply wraps around the R `SeuratDisk` package to perform the stated conversion. Included her for ease of use for the user.
#' Each assay in the Seurat Object is written into separate .h5ad files.
#'
#' @param SeuratObj Seurat Object
#' @param assays A vector of assay names from SeuratObj to be exported
#' @param h5ad_filename Filename of the output .h5ad file. Note that the final output filenames will have the assay names appended to this (see examples).
#'
#' @return A vector of .h5ad filenames which are outputted. Each file correspond to one Seurat assay, as indicated in the suffix inside the filename (see examples).
#'
#' @importFrom Seurat Assays
#' @importFrom SeuratDisk SaveH5Seurat Convert
#'
convertSeuratToH5ad <- function(SeuratObj, assays, h5ad_filename){
  # first convert all metadata columns which are factors into characters
  # this is to preserve them in the anndata object (if they are factors only the levels
  # will be retained, not the labels)
  for(col in colnames(SeuratObj@meta.data)){
    if(class(SeuratObj@meta.data[, col]) == "factor"){
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
    out_files <- c(out_files, out_file)
  }
  return( out_files )
}

#' parse the substring inside a given cell identifier which corresponds to the nucleotide barcode
#'
#' @description
#' `guessBarcodes` parses given cell identifiers to identify the substring which correspond to the nucleotide barcode included in the experiment.
#'
#' @details
#' Numeric / string prefices/suffices were typically added to cell identifiers to avoid wrong mapping across samples; however often these create issues when trying
#' to merge data on the *same* sample but annotated using different workflows.
#' This function attempts to resolve such issues by extracting the nucleotide barcodes actually introduced in the experiment.
#'
#' @param cell_name character, a cell identifier, typicall with prefix and/or suffix (e.g. "ACTGATGCAT-1", "SampleA_ATGAACCTATGG")
#' @param min_barcode_name minimum length of the nucleotide barcode (Default: 6)
#'
#' @return a vector with the input `cell_name` decomposed into these three entries:
#' \describe{
#'   \item{prefix} {prefix which exists in the input `cell_name` (NA if doesn't exist in `cell_name`)}
#'   \item{cell_name} {the actual nucleotide barcode}
#'   \item{suffix} {suffix which exists in the input `cell_name` (NA if doesn't exist in `cell_name`)}
#' }
#' @importFrom stringr str_detect str_locate str_sub
#'
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
#' `combineLoomFiles` combines .loom files generated using velocyto, on multiple BAM files,
#' into one loom file with the cell barcodes fixed to reflect the cell names in the given Seurat object.
#'
#' @details
#' `combineLoomFiles` take a vector of `sample_names` (which is assumed to be of the same length and
#' in the same order as `loom_files`), parse the prefix and suffix added to the cell barcodes belonging to the
#' given sample, and modify the column names of the matrices in the loom files accordingly. This allows for
#' combining the matrices coming from different samples without ambiguity of cell barcodes.
#' It then checks the matrices for overlap of genes with the given Seurat Object, and remove duplicated genes,
#' and finally write out the merged loom matrices into a new loom file, to be used for RNA velocity analysis.
#'
#' @param loom_files vector of loom files to be merged
#' @param new_loom_file, character, name of the new loom file to be written out containing the merged data
#' @param SeuratObj corresponding Seurat Object
#' @param sample_names vector of sample names to be looked for in the `seurat_sample_column` column of the Seurat metadata.
#' Assumed this is in the same order as the order given in `loom_files`.
#' @param seurat_sample_column column name in the Seurat metadata where the different values given in `sample_names` can be found
#'
#' @return a output message indicating success of writing out the merged loom matrices into the file given by `new_loom_file`.
#'
#' @importFrom velocyto.R read.loom.matrices
#' @import hdf5r
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
  loom_objs <- lapply(loom_files, velocyto.R::read.loom.matrices)
  # strip the -[0-9] suffix in the VDJ data frame barcode column
  barcodes_loom <- list()
  message("Assuming the order in sample_names correspond to the order in loom_files. If this is not the case please rerun this function ensuring the orderings of these vectors match up.")
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
#' `mergeVelocytoWithGEX` merges the velocyto spliced/unspliced gene counts with the AnnData object holding single-cell gene expression data.
#' This is the preprocessing function before calculating RNA velocity using the python scVelo package and workflow.
#'
#' @details
#' `mergeVelocytoWithGEX` uses the R reticulate package to run python commands which merge the velocyto counts with the anndata object.
#' It is assumed that the cell barcodes in the velocyto loom and the anndata object matches up. Please consult `guessBarcodes` and `combineLoomFiles` first to ensure these barcodes do match up.
#'
#' @param anndata_file filename pointing to the AnnData file containing gene expression data.
#' @param loom_file, character, name of the new loom file to be written out containing the merged data
#' @param anndata_out_filename output filename of the merged AnnData object to be written containing both the velocyto data and the gene expression data.
#' @param conda_env character, if not NULL this named conda environment is used to perform the merge.
#' (Default: NULL, i.e. no conda environment will be used, the program assumes the python packages `scanpy` and `scvelo` are installed in the local python)
#'
#' @return a output message indicating success of writing out the merged AnnData object into the file given by `anndata_out_filename`.
#'
#' @import reticulate
mergeVelocytoWithGEX <- function(anndata_file, loom_file, anndata_out_filename,
                                 conda_env = NULL)
{
  use_condaenv(conda_env)
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
#' `run_scVelo` calculates RNA velocity using the python scvelo package and the velocity models it implements.
#'
#' @details
#' `run_scVelo` uses the R reticulate package to run python commands which run scvelo RNA velocity calculations.
#' It follows the ["RNA Velocity Basics"](https://scvelo.readthedocs.io/VelocityBasics/) tutorial in the scvelo
#' documentation. Unfortunately due to conflicts of the plotting functionalities of R and python this function
#' does **NOT** implement the visualisation of velocity stream onto the dimensionality-reduced projection.
#'
#' @param anndata_file filename pointing to the AnnData file containing gene expression data and merged velocyto spliced/unspliced gene counts.
#' @param anndata_out_filename output filename of the merged AnnData object to be written the fitted RNA velocity estimates calculated using scvelo.
#' @param conda_env character, if not NULL this named conda environment is used to perform the merge.
#' (Default: NULL, i.e. no conda environment will be used, the program assumes the python packages `scanpy` and `scvelo` are installed in the local python)
#' @param scvelo_mode the 'mode' parameter in the python scvelo function `scv.tl.velocity`. (Default: "dynamical")
#' @param reduction the dimensionality reduction to project RNA velocity estimates onto (Default: "UMAP")
#' @param min_shared_counts include only genes detected in at least this number of cells. (Default: 20)
#' @param n_top_genes include only this many genes with the largest dispersion in the dataset (Default: 2000)
#'
#' @return a output message indicating success of writing out the AnnData object with merged scVelo results into the file given by `anndata_out_filename`.
#'
#' @import reticulate
run_scVelo <- function(anndata_file, anndata_out_filename, conda_env = NULL,
                       scvelo_mode = "dynamical", reduction = "UMAP",
                       min_shared_counts = 20, n_top_genes = 2000)
{
  use_condaenv(conda_env)
  py_run_string("import scanpy as sc")
  py_run_string("import scvelo as scv")
  py_run_string("scv.settings.verbosity = 3")
  py_run_string(paste0("adata = sc.read_h5ad('", anndata_file, "')"))
  # first calculate moments
  message("running scVelo pipeline ...")
  py_run_string(paste0("scv.pp.filter_and_normalize(adata, min_shared_counts=",
                       min_shared_counts, ", n_top_genes=", n_top_genes, ")"))
  py_run_string("scv.pp.neighbors(adata)")
  py_run_string("scv.pp.moments(adata, n_pcs=None, n_neighbors=None)")
  py_run_string("scv.tl.recover_dynamics(adata)")
  py_run_string(paste0("scv.tl.velocity(adata, mode = '", scvelo_mode, "')"))
  py_run_string("scv.tl.velocity_graph(adata)")
  py_run_string(paste0("scv.tl.velocity_embedding(adata, basis='",
                       stringr::str_to_lower(reduction), "', autoscale=False)"))
  py_run_string(paste0("adata.write_h5ad('", anndata_out_filename, "')"))
  return( paste0("AnnData object with scVelo results written to file '",
                 anndata_out_filename, "'.") )
}

#' Split AnnData object by levels in a specified meta data trait
#'
#' @description
#' `splitAnnData` splits a AnnData object by a given metadata column and write out the splitted data subsets into separate .h5ad files.
#'
#' @details
#' `splitAnnData` uses the R reticulate package to run python commands which subset the data by the metadata column
#' ( i.e. column in adata.obs) and write out each susbet as separate AnnData objects in .h5ad files.
#'
#' @param anndata_file filename pointing to the AnnData file.
#' @param split.by column name in the metadata (i.e. in `SeuratObj@meta.data` considering the Seurat object, or `adata.obs` considering the AnnData object)
#' indicating distinct levels to split the AnnData.
#' @param levels vector of values which can be found in the column given by `split.by`. Subsets of the Anndata by each value of `levels` will be made and written out into .h5ad files.
#' @param conda_env character, if not NULL this named conda environment is used to perform the merge.
#' (Default: NULL, i.e. no conda environment will be used, the program assumes the python packages `scanpy` and `scvelo` are installed in the local python)
#'
#' @return a vector of output .h5ad filenames with the indicated subset of the data (given in the filenames, see examples.)
#'
#' @import reticulate
splitAnnData <- function(anndata_file, split.by, levels, conda_env)
{
  # split the AnnData object by a given column ('split.by') in Anndata.obs,
  # given the set of values ('levels') to subset for.
  # Write out separate .h5ad file of the subsetted AnnData object
  use_condaenv(conda_env)
  py_run_string("import scanpy as sc")
  py_run_string(paste0("adata = sc.read_h5ad('", anndata_file, "')"))
  out_fn <- c()
  for(lvl in levels){
    message(paste0("Subsetting AnnData for ", split.by,  " == '", lvl, "' ..."))
    py_run_string(paste0("adata_subset = adata[adata.obs['", split.by, "'] == '", lvl, "']"))
    subset_filename <- paste0(gsub(".h5ad", "", basename(anndata_file)), "_", lvl, ".h5ad")
    subset_filename <- paste0(dirname(anndata_file), "/", fs::path_sanitize(subset_filename))
    py_run_string(paste0("adata_subset.write_h5ad('", subset_filename, "')"))
    out_fn <- c(out_fn, subset_filename)
  }
  return( out_fn )
}

#' Fit transition model on data using the python cellrank package
#'
#' @description
#' `fitTransitionModel` fits transition models on the data using the python cellrank package.
#'
#' @details
#' `fitTransitionModel` currently implements either the velocity kernel in cellrank (i.e. uses RNA velocity
#' information to fit transition probabilities) or the pseudotime kernel; a user-indicated column in the metadata
#' will be used as pseudotime reference to fit transition probabiltiies.
#' **NOTE:** In cases where the warning 'Biased KNN graph is disconnected' which will subsequently cause the `fitTPT` function in the pipeline to falied with error, in our experience it is likely to be caused by subsetting the data prior to computing transitions. Try setting `do_pca = FALSE` will preserve the original PCA and avoid this error.
#'
#' @param anndata_file filename pointing to the AnnData file.
#' @param conda_env character, if not NULL this named conda environment is used to perform the merge.
#' (Default: NULL, i.e. no conda environment will be used, the program assumes the python packages `scanpy`, `scvelo` and `cellrank` are installed in the local python)
#' @param mode character, either 'pseudotime' (uses the cellrank 'PseudotimeKernel') or 'velocity' (cellrank 'VelocityKernel'). (Default: 'pseudotime')
#' @param pseudotime_key character, column name which indicates the ranking to be used as pseudotime ordering of the cells. Not considered if mode is 'velocity'. (Default: 'csr_pot')
#' @param do_pca Should principal component analysis (PCA) be re-computed on the data? (Default: TRUE)
#'
#' @return a list with three entries:
#' \describe{
#'   \item{cellrank_obj}{Python `cellrank.tl.estimators.CFLARE` object containing details of the fitted transition model}
#'   \item{transition_matrix}{The matrix holding cell-to-cell transition probability (i.e. `cellrank_obj.transition_matrix`), but converted into a dense matrix. This is the single-cell transition matrix for downstream uses.}
#'   \item{CellID}{a vector of cell identifiers in the order of each row/column of `transition_matrix`.}
#' }
#'
#' @import reticulate
fitTransitionModel <- function(anndata_file, conda_env = NULL, mode = 'pseudotime',
                               pseudotime_key = 'csr_pot', do_pca = TRUE, do_neighbors = TRUE)
{
  if( ! mode %in% c("pseudotime", "velocity"))
    stop("Currently 'mode' must be either 'pseudotime' (cellrank PseudotimeKernel) or 'velocity' (cellrank VelocityKernel).")
  use_condaenv(conda_env)
  py_run_string("import scanpy as sc")
  py_run_string(paste0("adata = sc.read_h5ad('", anndata_file, "')"))
  if( do_pca ){
    py_run_string("sc.tl.pca(adata)")
  }
  if( do_neighbors ){
    py_run_string("sc.pp.neighbors(adata)")
  }
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
#' `fitTPT` fits Transition Path Theory on the transition model defined using `fitTransitionModel`, at
#' a 'coarse-grained' level where transitions are considered between *groups* of cells with grouping indicated by the user.
#'
#' @details
#' `fitTPT` interfaces with (and reimplements some routines to improve efficincy) the Python `deeptime` package to fit transition path theory (TPT) onto the markov state model defined by running
#' the `fitTransitionModel` function that uses cellrank under the hood. With the parameter `group.cells.by`, the user specifies
#' a scheme to group individual row/columns of the transition matrix (for example, by cell type or by isotype). The function
#' then fits TPT on to this grouped/'coarse-grained' transition matrix, upon user indicating a likely 'source' and 'target' state.
#' The output are estimated information flows ('flux') between different states in order to flow from the source to the target,
#' and the probabilities of sampling each state at equilibrium ('stationary distribution').
#' A random 'null background' model was fitted by randomly reshuffling columns of the transition matrix by `random_n` (default: 100) times.
#' These random fluxes help determine the significance of an observed flux, by calculating one-sided empirical probabilitys of the observed flux larger than the that observed in the randomised models.
#'
#' @param anndata_file filename pointing to the AnnData file.
#' @param CellrankObj the `cellrank_obj` entry of the output list of `fitTransitionModel`.
#' @param conda_env character, if not NULL this named conda environment is used to perform the merge.
#' (Default: NULL, i.e. no conda environment will be used, the program assumes the python packages `scanpy`, `scvelo` and `cellrank` are installed in the local python)
#' @param group.cells.by character, column in the metadata to group cells
#' @param source_state character, a value in the `group.cells.by` column which is taken as the source state for fitting transition path theory. All cells belonging to this group are considered as the source.
#' @param target_state character, a value in the `group.cells.by` column which is taken as the target state for fitting transition path theory. All cells belonging to this group are considered as the target.
#' @param random_n number of times to reshuffle transition matrix columns to derive randomised models (default: 100).
#'
#' @return a list with three entries:
#' \describe{
#'   \item{gross_flux}{a n-by-n matrix (where n is the total number of states), of total fluxes estimated between from a state (row) to another state (column). }
#'   \item{net_flux}{a n-by-n matrix (where n is the total number of states), of net fluxes estimated between from a state (row) to another state (column). (i.e. the difference between elements (i, j) and (j, i) of the gross flux matrix.) }
#'   \item{pathways}{a data.frame indicating the possible paths to take from `source_state` to `target_state`, and the likelihood (max: 100) to travel through each stated path.}
#'   \item{significance}{a n-by-n matrix (where n is the total number of states), where the observed gross flux is greater than the flux estimated in the randomised models. }
#'   \item{total_gross_flux}{element-wise sum of the gross_flux matrix. }
#'   \item{total_gross_flux_reshuffled}{element-wise sum of the gross_flux matrix, calculated over each randomised (randomly reshuffled transition matrix coluns) models.}
#'   \item{total_net_flux}{element-wise sum of the net_flux matrix. }
#'   \item{total_net_flux_reshuffled}{element-wise sum of the net_flux matrix, calculated over each randomised (randomly reshuffled transition matrix coluns) models.}
#'   \item{mfpt}{Mean First Passage Time required to travel from `source_state` to `target_state` as estimated by Transition Path Theory. }
#'   \item{mfpt_reshuffled}{Mean First Passage Time required to travel from `source_state` to `target_state` as estimated by Transition Path Theory, calculated over each randomised (randomly reshuffled transition matrix coluns) models. }
#'   \item{stationary_distribution}{Equilibrium probability of each state as estimated by Transition Path Theory. }
#'   \item{stationary_distribution_reshuffled}{Equilibrium probability of each state as estimated by Transition Path Theory, calculated over each randomised (randomly reshuffled transition matrix coluns) models. }
#' }
#'
#' @import reticulate
fitTPT <- function(anndata_file, CellrankObj, conda_env,
                   group.cells.by, source_state, target_state, random_n = 100)
{
  use_condaenv(conda_env)
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
  tpt[["pathways"]] <- reticulate::py_to_r(tpt[["pathways"]])
  tpt[["stationary_distribution_bootstrapping"]] <- lapply(tpt[["stationary_bootstrapping"]], unlist)
  tpt[["total_gross_flux"]] <- tpt[["total_gross_flux"]]
  tpt[["total_gross_flux_reshuffled"]] <- unlist(tpt[["total_gross_flux_randomised"]])
  tpt[["total_net_flux"]] <- tpt[["total_net_flux"]]
  tpt[["total_net_flux_reshuffled"]] <- unlist(tpt[["total_net_flux_randomised"]])
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
  dimnames(tpt[["net_flux"]]) <- list(cluster_order, cluster_order)
  tpt[["net_flux"]] <- tpt[["net_flux"]][cluster_names, cluster_names]
  dimnames(tpt[["significance"]]) <- list(cluster_order, cluster_order)
  tpt[["significance"]] <- tpt[["significance"]][cluster_names, cluster_names]
  return(tpt[c("gross_flux", "net_flux", "pathways", "significance",
               "total_gross_flux", "total_gross_flux_reshuffled",
               "total_net_flux", "total_net_flux_reshuffled",
               "mfpt", "stationary_distribution",
i              "mfpt_reshuffled", "stationary_distribution_reshuffled", 
               "stationary_distribution_bootstrapping")])
}

#' Find the centroid of each cell cluster in the dimensionality-reduced space.
#'
#' @description
#' `FindClusterCentroid` calculates for each cell cluster, the centroid coordinates in a specified dimensionality reduction space.
#'
#' @details
#' `FindClusterCentroid` calculates the mean coordinate of the 2 axes of the dimensionality reduction space
#' specified by `dim_reduce`, for each cluster specified by `group.by`.
#'
#' @param SeuratObj Seurat object
#' @param group.by column in Seurat object meta.data to group the cells (Default: 'seurat_clusters')
#' @param dim_reduce name of the dimensionality reduction to calculate centroid coordinates (Default: "UMAP")
#'
#' @return a data.frame of the coordinates (in the 2 axes of the dimenesionality reduced space) for each cell cluster.
#'
#' @importFrom stringr str_to_upper
#' @importFrom plyr ddply
#' @importFrom Seurat Reductions
FindClusterCentroid <- function(SeuratObj, group.by = "seurat_clusters", dim_reduce = "UMAP")
{
  dim_reduce <- stringr::str_to_upper(dim_reduce)
  if( ! dim_reduce %in% Seurat::Reductions(SeuratObj) )
    stop("'dim_reduce' is not found in the list of dimensionality reductions available for the given Seurat object.")
  # calculate cluster centroid (i.e. mean UMAP_1, UMAP_2)
  d <- Seurat::FetchData(SeuratObj, vars = c(paste0(dim_reduce, "_1"), paste0(dim_reduce, "_2"), group.by))
  o <- plyr::ddply(d, .variables = group.by, function(x){
    c(paste0(dim_reduce, "_1") = mean(x[, paste0(dim_reduce, "_1")]),
      paste0(dim_reduce, "_2") = mean(x[, paste0(dim_reduce, "_2")]))
  })
  oo <- as.matrix(o[, c(paste0(dim_reduce, "_1"), paste0(dim_reduce, "_2"))])
  rownames(oo) <- o[, group.by]
  return(oo)
}

#' Visualise fluxes between clusters with a network plot
#'
#' @description
#' `plotFlux` plots fluxes between cell clusters in the form of a network plot.
#'
#' @details
#' `plotFlux` represents fluxes between cell clusters in a graph, where the thickness of the edges signifies the
#' amount of fluxes between the given pair of clusters. The placement of clusters respect their placements in the
#' default dimensionality-reduced projection of the gene expression data stored in the Seurat Object (e.g. the UMAP
#' projection calculated using Seurat). It can either be a stand-alone plot, or be added on top of UMAP projections (see example).
#'
#' @param TPTObj List of TPT results, output from the `fitTPT` function.
#' @param SeuratObj Seurat object
#' @param significance_threshold The minimum value in the `TPTObj$significance` matrix for a flux to be shown in the
#' resulting plot (all fluxes with significance below this value will be removed from the visualisation).
#' @param mask_lower_tri Should fluxes which are considered 'reversed' (i.e. in the lower triangle of the flux matrix)
#' be masked? (Default: FALSE, i.e. all fluxes in the gross_flux matrix is shown)
#' @param mask_threshold the minimum percentage (max 100) of total flux to be shown in the resulting plot. (Default: 1, i.e.
#' gross fluxes below 1 will be removed from the visualisation)
#' @param new_plot Should a stand-alone plot be returned. If FALSE, ggplot2 geoms are returned for adding on top of e.g. Seurat::DimPlot. (Default: FALSE)
#'
#' @return If `new_plot` is TRUE, a network plot where the nodes are placed respecting the positioning of cell clusters in the dimensionality-reduction
#' projection of the data, and the thickness of edges connecting the nodes signify the flux between pairs of clusters. If `new_plot` is FALSE,
#' the `ggplot2` geom objects are returned to enable plotting the network graph on top of dimensionality projections (see example, where the UMAP
#' projection is plotted using the function `dim_plot` implemented in this package).
#'
#' @import ggplot2
#' @importFrom stringr str_to_upper
#' @importFrom plyr ddply
#' @importFrom Seurat Reductions
plotFlux <- function(TPTObj, SeuratObj,
                     significance_threshold = 0.9,
                     mask_lower_tri = FALSE, mask_threshold = 1,
                     new_plot = TRUE)
{
  flux_matrix <- TPTObj$gross_flux
  significance_matrix <- TPTObj$significance
  stationary_distribution <- TPTObj$stationary_distribution
  source("/media/josephn/Seagate4TB/GLT_datasets/ggplot2_multiple_scales.R")
  library(ggnetwork)
  flux_matrix[significance_matrix < significance_threshold] <- 0
  if(mask_lower_tri){
    flux_matrix[lower.tri(flux_matrix)] <- 0
  } else if( is.numeric(mask_threshold) ) {
    flux_matrix[flux_matrix <= mask_threshold] <- 0
  }
  graph <- igraph::graph_from_adjacency_matrix(flux_matrix,
                                               weighted = TRUE, mode = "directed")
  graph <- igraph::set_vertex_attr(graph, "nodeweight", value = stationary_distribution)
  # initialise ggnetcombineLoomFileswork
  cluster_pos <- FindClusterCentroid(SeuratObj)
  graph <- ggnetwork(graph, layout = cluster_pos[rownames(flux_matrix), ])
  if( new_plot ){
    ggplot(graph,
           aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(color = "grey50",
                 arrow = arrow(length = unit(6, "pt"), type = "closed"),
                 curvature = 0.2, aes(size = weight)) +
      scale_size_continuous(range = c(0.5, 2)) +
      new_scale("size") +
      geom_point(aes(size = nodeweight * 2), color = 'orange') +
      geom_nodetext(aes(label = name, size = nodeweight)) +
      scale_size_continuous(range = c(3, 8)) +
      theme_blank() + theme(legend.position = "none")
  } else {
    # rescale x, y, xend, yend
    graph[, "x"] <- graph[, "x"] * abs(diff(range(cluster_pos[, "UMAP_1"]))) + min(cluster_pos[, "UMAP_1"])
    graph[, "xend"] <- graph[, "xend"] * abs(diff(range(cluster_pos[, "UMAP_1"]))) + min(cluster_pos[, "UMAP_1"])
    graph[, "y"] <- graph[, "y"] * abs(diff(range(cluster_pos[, "UMAP_2"]))) + min(cluster_pos[, "UMAP_2"])
    graph[, "yend"] <- graph[, "yend"] * abs(diff(range(cluster_pos[, "UMAP_2"]))) + min(cluster_pos[, "UMAP_2"])
    list(
      geom_edges(color = "grey50",
                 arrow = arrow(length = unit(6, "pt"), type = "closed"),
                 curvature = 0.2, data = graph,
                 aes(x = x, y = y, xend = xend, yend = yend, size = weight)),
      scale_size_continuous(range = c(0.5, 2), guide = "none"),
      new_scale("size"),
      geom_point(data = graph, aes(x = x, y = y, size = nodeweight * 2),
                 color = 'orange'),
      geom_nodetext(data = graph, aes(x = x, y = y, label = name)),
      scale_size_continuous(range = c(0.1, 2), guide = "none")

    )
  }
}

#' Plot dimensionality-reduced projection of single-cell data
#'
#' @description
#' `dim_plot` plots the dimensionality-reduced projection stored in the Seurat object.
#'
#' @details
#' `dim_plot` effectively emulates the `DimPlot` function in Seurat, but here returns a native `ggplot2` object
#' whch enables additional geoms (e.g. the flux network plot from `plotFlux`) to be drawn on top of the dimensionality-reduced
#' projection.
#'
#' @param SeuratObj Seurat object
#' @param group.by column in Seurat object meta.data to group the cells (Default: 'seurat_clusters')
#' @param dim_reduce name of the dimensionality reduction to calculate centroid coordinates (Default: "UMAP")
#'
#' @return ggplot2 object of the dimensionality-reduced projection of the single-cell data, on the first 2 axes
#' of the projection given in `dim_reduce`
#'
#' @import ggplot2
#' @importFrom stringr str_to_upper
dim_plot <- function(SeuratObj, group.by = "seurat_clusters", dim_reduce = "UMAP")
{
  # emulate Seurat::DimPlot
  dim_reduce <- stringr::str_to_upper(dim_reduce)
  d <- Seurat::FetchData(SeuratObj, c(paste0(dim_reduce, "_1"), paste0(dim_reduce, "_2"), group.by))
  n_clusters <- length(levels(d[, group.by]))
  ggplot(d, aes_string(x = paste0(dim_reduce, "_1"), y = paste0(dim_reduce, "_2"))) +
    geom_point(aes_string(colour = group.by)) +
    scale_color_manual(values = scales::hue_pal()(n_clusters),
                       name = group.by) +
    cowplot::theme_cowplot() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(), axis.title = element_blank())
}

#' Prepare class-switch recombination (CSR) transition data for visualisation
#'
#' @description
#' `prepareCSRtransitions` parse data from `fitTPT` to prepare them for visualising the estimated isotype-switching dynamics.
#'
#' @details
#' `prepareCSRtransitions` performs the following processing steps to parse data from the `fitTPT` function for visualising them using `plotCSRtransitions`:
#' \itemize{
#'   \item{filter flux matrix by significance}{the gross_flux matrix from Transition Path Theory (TPT) is filtered using the significance calculated empirically upon
#'   comparison with fluxes generated in randomised models. Only fluxes greater than the randomised fluxes (threshold given by `significance_threshold`, default 0.9) are shown.}
#'   \item{filter flux matrix by magnitude}{gross_flux matrix is filtered further such that only those fluxes representing more than a given percentage of total fluxes are shown. Default is 1%.}
#'   \item{remove improbable CSR combinations}{isotype switches which are improbable (i.e. switching back to an isotype 5' to the current) are removed. This can be turned off by indicating `mask_improbable_csr = FALSE`.}
#' }
#' A parsed output is returned to be directly fed into `plotCSRtransitions` to visualise the dynamics of CSR in the given dataset estimated using productive/sterile transcript information.
#'
#' @param TPTObj List of TPT results, output from the `fitTPT` function.
#' @param SeuratObj Seurat object
#' @param ighc_count_assay_name name of assay in SeuratObj which holds the IgH productive/sterile transcript counts. (Default: "IGHC")
#' @param significance_threshold The minimum value in the `TPTObj$significance` matrix for a flux to be shown in the
#' resulting plot (all fluxes with significance below this value will be removed from the visualisation).
#' @param mask_improbable_csr Should isotype combinations which represents improbable Class-switch recombination events (i.e. switching back to an isotype 5' to the current one) be removed from visualisation? (Default: TRUE)
#' @param mask_threshold the minimum percentage (max 100) of total flux to be shown in the resulting plot. (Default: 1, i.e.
#' gross fluxes below 1 will be removed from the visualisation)
#'
#' @return A named list of three entries:
#' \describe{
#'   \item{stationary_distribution}{A data.frame of isotypes, their x/y position in the plot, and their stationary distribution estimated from TPT.}
#'   \item{flux}{A data.frame listing the fluxes which survive the filtering steps applied (see 'Details'). Their x/y positions in the plot (`from`/`to`) and their gross fluxes are listed.}
#'   \item{c_genes}{a vector of gene names representing the IgH C genes considered in the model.}
#' }
#' This list can be directly passed to the `plotCSRtransitions` function for visualisation.
#'
#' @importFrom reshape2 melt
#' @importFrom stringr str_to_lower str_replace str_to_sentence
prepareCSRtransitions <- function(TPTObj, SeuratObj,
                                  ighc_count_assay_name = "IGHC",
                                  significance_threshold = 0.9,
                                  mask_improbable_csr = TRUE, mask_threshold = 1)
{
  library(ggplot2)
  flux_matrix <- TPTObj$gross_flux
  significance_matrix <- TPTObj$significance
  stationary_distribution <- TPTObj$stationary_distribution
  bs_stationary <- TPTObj$stationary_distribution_bootstrapping
  total_flux <- TPTObj$total_gross_flux

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
      stationary_distribution <- c(stationary_distribution, 0)
      names(stationary_distribution)[length(stationary_distribution)] <- c_gene
    }
  }
  flux_matrix <- flux_matrix[c_genes, c_genes]
  significance_matrix <- significance_matrix[c_genes, c_genes]
  # filtering
  flux_matrix[significance_matrix < significance_threshold] <- 0
  if(mask_improbable_csr){
    flux_matrix[lower.tri(flux_matrix)] <- 0
  }
  if( is.numeric(mask_threshold) ) {
    flux_matrix[flux_matrix <= mask_threshold] <- 0
  }
  graph <- reshape2::melt(flux_matrix, value.name = "flux", varnames = c("from", "to"))
  graph$from <- as.numeric(graph$from)
  graph$to <- as.numeric(graph$to)
  graph$flux <- graph$flux * total_flux / 100
  graph <- graph[which(graph$flux > 0), ]
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
  # weight the fluxes by stationary distribution
  # graph <- merge(graph, isotypes[, c("x", "stationary")], by.x = "from", by.y = "x",
  #                all.x = TRUE, all.y= FALSE, sort = FALSE)
  # colnames(graph)[ncol(graph)] <- "stationary_from"
  # graph <- merge(graph, isotypes[, c("x", "stationary")], by.x = "to", by.y = "x",
  #                all.x = TRUE, all.y= FALSE, sort = FALSE)
  # colnames(graph)[ncol(graph)] <- "stationary_to"
  # graph$flux <- graph$flux * graph$stationary_from * graph$stationary_to
  return(list("stationary_distribution" = isotypes, "flux" = graph[, c("from", "to", "flux")],
              "c_genes" = c_genes))
}

#' Visualising class-switch recombination (CSR) transitions estimated from the data
#'
#' @description
#' `plotCSRtransitions` visualises the isotype-switching dynamics estimated on the given data using transition models implemented in this package.
#'
#' @details
#' `plotCSRtransitions` either takes the output of `prepareCSRtransitions` (argument `csr_transitions`) and generate a
#' bar plot showing stationary distribution of the isotypes and arrows indicating the amount of fluxes estimated from
#' the fitted transition model. Alternatively, user can supply the fitted Transition Path Theory object (output of the
#' `fitTPT` function) as argument `TPTObj`, and the associated `SeuratObj`; in this case the `prepareCSRtransitions` function
#' will be called internally to parse the transition data for visualisation.
#'
#' @param csr_transitions list of parsed transitions from the `prepareCSRtransitions` package. If this is supplied and not NULL, the `TPTObj` and `SeuratObj` arguments will be ignored.
#' @param TPTObj List of TPT results, output from the `fitTPT` function. Considered only if `csr_transitions` is NULL.
#' @param SeuratObj Seurat object. Considered only if `csr_transitions` is NULL.
#' @param ighc_count_assay_name name of assay in SeuratObj which holds the IgH productive/sterile transcript counts. (Default: "IGHC")
#' @param return_plot Should the CSR transition plot be returned? If FALSE, a named list of `stationary_distribution` and `flux` will be returned which contains the data frames to be visualised in this plot. (Default: TRUE)
#' @param significance_threshold The minimum value in the `TPTObj$significance` matrix for a flux to be shown in the
#' resulting plot (all fluxes with significance below this value will be removed from the visualisation).
#' @param mask_improbable_csr Should isotype combinations which represents improbable Class-switch recombination events (i.e. switching back to an isotype 5' to the current one) be removed from visualisation? (Default: TRUE)
#' @param mask_threshold the minimum percentage (max 100) of total flux to be shown in the resulting plot. (Default: 1, i.e.
#' gross fluxes below 1 will be removed from the visualisation)
#' @param curvature amount of curvature of the arrows representing CSR fluxes. (Default: 0.1)
#' @param arrow `grid::arrow` object specifying the size and aesthetics of the plotted arrows representing CSR fluxes.
#' @param arrow_colour Optional. a column from the 'flux' data frame holding numeric data to be visualised as heat scale in the arrow colours. Useful for
#' overlaying custom comparisons (e.g. flux differences between Wild-type and Knock-down conditions) onto the plot. If NULL, use the flux amounts to scale
#' the opaqueness of the black arrows. (Default: NULL)
#' @param bar_colour optiona. a column from the 'stationary_distribution' data frame holding categorical information to show stationary distribution from differnet conditions as different coloured bars in the bar plot.
#' If NULL, only 1 bar will be shown per isotype, in grey. (Default: NULL)
#'
#' @return A ggplot2 object containing a bar plot showing stationary distribution of the isotypes, and arrows
#' indicating the amount of class-switch recombination events estimated from the fitted transition model
plotCSRtransitions <- function(csr_transitions = NULL,
                               TPTObj = NULL,
                               SeuratObj = NULL,
                               ighc_count_assay_name = "IGHC",
                               return.plot = TRUE,
                               significance_threshold = 0.9,
                               mask_improbable_csr = TRUE, mask_threshold = 1,
                               curvature = 0.1,
                               arrow = grid::arrow(length = unit(6, "pt"),
                                                   type = "closed"),
                               arrow_colour = NULL, bar_colour = NULL)
{
  if( !is.null(csr_transitions) ){
    if( sum(is.na(match(names(csr_transitions), c("stationary_distribution", "flux", "c_genes")))) > 0 ){
      stop("'csr_transitions' appear not to be derived from the output of 'prepareCSRtransitions'.\n
           If you start from a TPT object and a Seurat object, indicate the argument names 'TPTObj' and 'SeuratObj' when you pass these objects to this plotCSRtransitions function.")
    }
    return(
      plotCSRtransitions_(csr_transitions, return.plot = return.plot,
                          curvature = curvature,
                          arrow = arrow, arrow_colour = arrow_colour,
                          bar_colour = bar_colour)
    )
  }
  prepared <- prepareCSRtransitions(
    TPTObj, SeuratObj,
    ighc_count_assay_name = ighc_count_assay_name,
    significance_threshold = significance_threshold,
    mask_improbable_csr = mask_improbable_csr, mask_threshold = mask_threshold
  )
  plotCSRtransitions_(prepared, return.plot = return.plot,
                      curvature = curvature,
                      arrow = arrow, arrow_colour = arrow_colour,
                      bar_colour = bar_colour)
}

#' Workhorse function for visualising class-switch recombination (CSR) transitions
#'
#' @description
#' `plotCSRtransitions_` is the workhorse for plotting the CSR transition data. Called
#' by `plotCSRtransitions` to generate the ggplot2 object.
#'
#' @details
#' Called by `plotCSRtransitions` to generate the ggplot2 object.
#'
#' @param prepared list of parsed transitions from the `prepareCSRtransitions` package. If this is supplied and not NULL, the `TPTObj` and `SeuratObj` arguments will be ignored.
#' @param return_plot Should the CSR transition plot be returned? If FALSE, a named list of `stationary_distribution` and `flux` will be returned which contains the data frames to be visualised in this plot. (Default: TRUE)
#' @param significance_threshold The minimum value in the `TPTObj$significance` matrix for a flux to be shown in the
#' resulting plot (all fluxes with significance below this value will be removed from the visualisation).
#' @param mask_improbable_csr Should isotype combinations which represents improbable Class-switch recombination events (i.e. switching back to an isotype 5' to the current one) be removed from visualisation? (Default: TRUE)
#' @param mask_threshold the minimum percentage (max 100) of total flux to be shown in the resulting plot. (Default: 1, i.e.
#' gross fluxes below 1 will be removed from the visualisation)
#' @param curvature amount of curvature of the arrows representing CSR fluxes. (Default: 0.1)
#' @param arrow `grid::arrow` object specifying the size and aesthetics of the plotted arrows representing CSR fluxes.
#' @param arrow_colour Optional. a column from the 'flux' data frame holding numeric data to be visualised as heat scale in the arrow colours. Useful for
#' overlaying custom comparisons (e.g. flux differences between Wild-type and Knock-down conditions) onto the plot. If NULL, use the flux amounts to scale
#' the opaqueness of the black arrows. (Default: NULL)
#' @param bar_colour optiona. a column from the 'stationary_distribution' data frame holding categorical information to show stationary distribution from differnet conditions as different coloured bars in the bar plot.
#' If NULL, only 1 bar will be shown per isotype, in grey. (Default: NULL)
#'
#' @return A ggplot2 object containing a bar plot showing stationary distribution of the isotypes, and arrows
#' indicating the amount of class-switch recombination events estimated from the fitted transition model. If
#' `return.plot == FALSE`, a list of 2 items containng the 'stationary_distribution' and 'flux' data frames
#' from `prepared`.
plotCSRtransitions_ <- function(prepared,
                                return.plot = TRUE,
                                curvature = 0.1,
                                arrow = grid::arrow(length = unit(6, "pt"), type = "closed"),
                                arrow_colour = arrow_colour, bar_colour = bar_colour)
{
  isotypes <- prepared[["stationary_distribution"]]
  graph <- prepared[["flux"]]
  c_genes <- prepared[["c_genes"]]
  params <- list(arrow = arrow, curvature = -curvature,
                 angle = 90, ncp = 5)
  arrow_pos <- isotypes[which(isotypes$isotype %in% c_genes[1:max(c(graph$from, graph$to))]),
                        "stationary"]
  arrow_pos <- max(arrow_pos)
  if( !is.null(bar_colour) ){
    p <- ggplot(isotypes) +
      geom_bar(aes_string(x = "isotype", y = "stationary", fill = bar_colour), 
               stat = "identity", position = position_dodge2()) +
      geom_errorbar(aes_string(x = "isotype", ymin = "lowq", ymax = "highq"), 
		    width = 0, position = position_dodge2())
  } else {
    p <- ggplot(isotypes) +
      geom_bar(aes_string(x = "isotype", y = "stationary"), 
               stat = "identity", position = position_dodge2()) +
      geom_errorbar(aes_string(x = "isotype", ymin = "lowq", ymax = "highq"),
                    width = 0, position = position_dodge2())
  }
  if( !is.null(arrow_colour) ){
    alpha <- abs(graph[, arrow_colour])/max(abs(graph[, arrow_colour]))
    p <- p + ggplot2::layer(data = graph,
                            mapping = aes_string(x = "from", xend = "to", y = arrow_pos,
                                                 yend = arrow_pos, color = arrow_colour,
                                                 alpha = alpha),
                            stat = "identity", position = position_nudge(y = 0.1),
                            geom = ggplot2::GeomCurve, params = params)
  } else {
    alpha <- abs(graph[, "flux"])/max(abs(graph[, "flux"]))
    p <- p + ggplot2::layer(data = graph,
                            mapping = aes_string(x = "from", xend = "to", y = arrow_pos,
                                                 yend = arrow_pos, alpha = alpha),
                            stat = "identity", position = position_nudge(y = 0.1),
                            geom = ggplot2::GeomCurve, params = params)
  }
  p <- p + scale_y_continuous(limits = c(0, 1), name = "stationary distribution") +
    scale_x_discrete(drop=FALSE) +
    geom_text(aes(x = isotype, y = max(stationary), label = isotype), vjust = -1) +
    scale_color_gradient2(name = arrow_colour) + scale_alpha_continuous(limits = c(0, 1), guide = "none") +
    cowplot::theme_cowplot() + #scale_x_discrete(breaks = rev(c_genes)) +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), axis.title.x = element_blank())
  if( arrow_pos > 0.7 ){
    suppressMessages(
      p <- (p + scale_y_continuous(limits = c(0, 1.3), name = "stationary distribution"))
    )
  }
  if( return.plot ) return(p)
  else return(list("stationary_distribution" = isotypes, "flux" = graph))
}

#' Computing distances between multiple transition matrices
#'
#' @description
#' `compareTransitionMatrices` computes distances between a list of transition matrices.
#'
#' @details
#' `compareTransitionMatrices` takes the list of transition matrices given in `matrix_list`, and sample realisations from
#' the Markov chain defined using the transition matrix. It then compares the similarity of these sampled trajectories
#' using either the Kullback-Leibler divergence (`distance_metric == "KL"`) or the Jensen-Shannon divergence (`distance_metric == "JSD"`).
#' Both can be interpreted the same way (the larger this number, the more different two transition matrices and their
#' realisations are), although the Jensen-Shannon divergence is scaled between 0 and 1.
#'
#' @param matrix_list list of transition matrices to be compared.
#' @param SeuratObj Seurat object. Considered only if `csr_transitions` is NULL.
#' @param cells vector of cell identifiers corresponding to the row/column order in the transition matrices. (Assumed that all supplied transition matrices have exactly the same row/column ordering.)
#' @param group.by column in the Seurat object metadata on which cells are grouped. (Default: 'seurat_clusters')
#' @param n_realisations number of trajectories ('realisations') to be sampled from the Markov model defined using each transition matrix. (Default: 1000)
#' @param n_step number of time-steps in each trajectory/realisation to be sampled. (Default: 1000)
#' @param distance_metric the distance metric to be calculated. Either "KL" (for Kullback-Leibler divergence) or "JSD" (Jensen-Shannon divergence). (Default: "KL")
#'
#' @return A list of two entries:
#' \describe{
#'   \item{distance}{distance (`distance_metric`) between the trajectories sampled from each possible pair of transition matrix.}
#'   \item{sampled_transitions}{list of trajectories sampled from the transition matrices.}
#' }
#'
compareTransitionMatrices <- function(matrix_list, SeuratObj,
                                      cells, group.by = "seurat_clusters",
                                      n_realisation = 1000, n_step = 1000,
                                      distance_metric = "KL")
{
  if( ! distance_metric %in% c("KL", "JSD") )
    stop("'distance_metric' must be one of 'KL' (Kullback-Leibler divergence, default) or 'JSD' (Jensen-Shannon divergence).")
  library(markovchain)
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
  dimnames(dist_mat) <- list(names(sampled_transitions), names(sampled_transitions))
  return(list("distance" = dist_mat,
              "sampled_transitions" = sampled_transitions))
}
