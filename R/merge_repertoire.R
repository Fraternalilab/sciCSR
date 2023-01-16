#' Annotate heavy-light chain pairing
#'
#' @description
#' `annotatePairing` checks the repertoire table to annotate heavy-light chain pairing. This function classifies for a given cell (i.e. a given cell barcode)
#' whether it is a singlet/doublet etc. based on the number of heavy (H) and light (L) chain sequences observed in the data
#'
#' @details
#' The function assumes input of sequences all belonging to the same cell (i.e. bearing the same cell barcode).
#' Based on the number of heavy and light chain sequences associated with a given cell barcode they are classified into one of the following categories:
#' \describe{
#'   \item{BCR_doublets}{More than 1 H and More than 1 L per cell barcode}
#'   \item{multi_LC_diff_class}{Exactly 1 H and more than 1 L per cell barcode, both kappa and lambda light chains are found}
#'   \item{multi_LC_same_class}{Exactly 1 H and more than 1 L per cell barcode, light chains are either all kappa or all lambda}
#'   \item{multi_HC}{More than 1 H and exacly 1 L per cell barcode}
#'   \item{singlet}{Exactly 1 H and/or exactly 1 L per cell barcode}
#' }
#'
#' @param tb_list list of data.frame holding VDJ data. Each element corresponds to subset of sequences bearing one specific cell barcode
#' @param c_gene_column Column name indicating where the isotype/light chain type information is stored in the data.frames in `tb_list`.
#'
#' @return Any one of 'singlet', 'BCR_doublet', 'multi_LC_same_class', 'multi_LC_diff_class' or 'multi_HC', according to the rules stated in the 'Details' section.
#' @export annotatePairing
annotatePairing <- function(tb_list, c_gene_column = 'c_gene')
{
  n_H <- nrow( tb_list[['H']] )
  n_L <- nrow( tb_list[['L']] )
  l_class <- tb_list[['L']][, c_gene_column]
  l_class <- unique(substr(l_class, 1, 4))
  if( 'Multi' %in% tb_list[['L']]$chain ) return('BCR_doublet')
  if( n_H > 1 && n_L > 1) return('BCR_doublet')
  if( n_H > 1 && n_L == 1) return('multi_HC')
  if( n_H == 1 && n_L > 1 && length(l_class) == 1)
    return('multi_LC_same_class')
  if( n_H == 1 && n_L > 1 && length(l_class) > 1)
    return('multi_LC_diff_class')
  return('singlet')
}

#' Collapse the VDJ repertoire data frames by cell barcode
#'
#' @description
#' `collapseBCR` considers an input data frame of scBCR-seq contig annotations and collapses this data frame such that for each cell at most 1 heavy chain and 1 ligh chain are retained as representative, in order to merge these data into the Seurat object holding scRNA-seq data where data are organised at a per-cell level.
#'
#' @details
#' The function cleans the barcode table such that potential doublets and/or cells with more than 1 heavy/light chain sequences are dealt with prior to merging with scRNA-seq data in the form of Seurat objects to avoid problems in merging owing to one-to-many mapping between cell barcode and sequences.
#' The most frequently observed contig is taken as the representative.
#'
#' @param tb a single data.frame holding VDJ data of 1 library
#' @param format data format. Default = '10X'. For the moment only '10X' is supported.
#' @param full.table if TRUE, will return a list with two elements (1) the collapsed data with max 1H and 1L per cell, and (2) the full table (Default: FALSE)
#'
#' @return If `full.table` is FALSE, an input data frame where for each cell at most 1 heavy chain and 1 light chain will be retained as representative. Otherwise, it will be a list of 2 data frames, one the collapsed version and one the full version.
#' An additional column 'bcr_type' will be added to the data frame indicating whether each cell barcode is a singlet/doublet etc. based on the number of observed heavy/light chains. (see `?annotatePairing()`)
#' @export collapseBCR

collapseBCR <- function(tb, format = '10X', full.table = FALSE)
{
  if( format == '10X' ){
    cell_id <- 'barcode'
    locus <- 'chain'
    umi_count <- 'umis'
  } else {
    stop("Only '10X' is allowed in the argument 'format'.")
  }
  tb <- split(tb, f = tb[, cell_id])
  tb <- lapply(tb, function(tb_part){
    h_or_l <- sapply(tb_part[, locus], function(x){
      if(x == "IGH") return('H')
      else if(x == 'IGL') return('L')
      else if(x %in% c('M', 'G', 'A', 'E', 'D')) return('H')
      return('L')
    })
    h_or_l <- factor(h_or_l, levels = c('H', 'L'))
    tb_part <- split(tb_part, f = h_or_l)
    names(tb_part) <- levels(h_or_l)
    bcr_type <- annotatePairing(tb_part)
    tb_part_collapsed <- do.call("rbind", lapply(tb_part, function(tb_subset){
      tb_subset[which.max(tb_subset[, umi_count]), ]
    }))
    tb_part_collapsed$bcr_type <- bcr_type
    if( full.table ){
      tb_part <- do.call("rbind", tb_part)
      tb_part$bcr_type <- bcr_type
      return( list('collapsed.table' = tb_part_collapsed, 'full.table' = tb_part) )
    } else return( tb_part_collapsed )
   })
  if( full.table ){
    return( list('collapsed.table' = do.call("rbind", lapply(tb, function(x) x[[1]])),
                 'full.table' = do.call("rbind", lapply(tb, function(x) x[[2]]))))
  } else {
    return( do.call( "rbind", tb ))
  }
}

#' Repairing cell barcodes in VDJ data frames to match the Seurat object
#'
#' @description
#' `repairBarcode` map samples to their specific barcode suffixes in the Seurat object and alters the cell barcode given in the VDJ data frame to match the Seurat object. This facilitates merging the repertoire data into the Seurat object as additional metadata columns.
#'
#' @details
#' The function first establishes a mapping table between samples (given in `seurat_sample_column` in the Seurat object) and their sample-specific cell barcode suffix, and then repairs the barcodes given in the VDJ data frame to match the Seurat object.
#'
#' @param vdj_list a list of data.frame, each holding VDJ data of 1 library
#' @param SeuratObj Seurat object.
#' @param sample_names a vector of sample names in the data. 'Sample' here refers to sequencing library.
#' @param seurat_sample_column Column name in `SeuratObj` metadata where the sample name/ID information is stored.
#'
#' @return A list of vdj data.frames with the cell barcodes in the column 'barcode' repaired to match barcodes in the Seurat object.
#' @export repairBarcode

repairBarcode <- function(vdj_list, SeuratObj, sample_names,
                          seurat_sample_column = 'sample_id')
{
  # Here establish a mapping table from the Seurat object and then
  # overwrite the barcodes in the VDJ data frames

  # strip the -[0-9] suffix in the VDJ data frame barcode column
  for (i in seq_along(vdj_list)) {
    vdj_list[[i]]$barcode <- sapply(vdj_list[[i]]$barcode, function(x){
      unlist(strsplit(x, split = "-"))[1]
    })
  }

  # mapping from suffix to sample names/IDs
  sample_map <- SeuratObj@meta.data[, c(seurat_sample_column, "orig.ident")]
  sample_map$suffix <- sapply(rownames(sample_map), function(x){
    unlist(strsplit(x, split = "-"))[2]
  })
  sample_map <- unique(sample_map[, c(seurat_sample_column, "suffix")])
  rownames(sample_map) <- sample_map[, seurat_sample_column]
  sample_map <- sample_map[sample_names, ]
  sample_map <- sample_map$suffix
  names( sample_map ) <- sample_names
  # overwrite cell barcodes by replacing suffixes
  for (i in seq_along(vdj_list)) {
    sample_id <- names(vdj_list)[i]
    vdj_list[[i]]$barcode <- paste0(vdj_list[[i]]$barcode,
                                    "-", sample_map[sample_id])
    vdj_list[[i]]$sample_name <- sample_id
  }

  vdj_list
}

#' Add annotations from VDJ data frame into the Seurat object
#'
#' @description
#' `combineBCR` merges annotations from data frames holding repertoire data as metadata of a Seurat object.
#'
#' @details
#' The function first reorganises the VDJ data frame such that one line represent 1 cell.
#' It then selects columns from this data frame and add them as metadata in the given Seurat object.
#'
#' @param vdj a data.frame of VDJ data after running `collapseBCR`.
#' @param SeuratObj Seurat object.
#' @param keep_columns a vector of column names in `vdj` to be kept and added to the metadata slot of the Seurat object.
#' By default the following informations are included: VDJ gene usage, isotype, CDR3 sequence, as well as binary indications of full-length/productive.
#' @param seurat_sample_column (Optional) Column name in `SeuratObj` metadata where the sample name/ID information is stored. Used when the barcode column from vdj data frame does not match cell names from the Seurat object given by `Cells(SeuratObj)`.
#' @param seurat_cell_name_column (Optional) Column name in `SeuratObj` metadata where the cell barcodes are stored. Used when the barcode column from vdj data frame does not match cell names from the Seurat object given by `Cells(SeuratObj)`.
#'
#' @return Seurat object with VDJ annotations added to the metadata slot.
#' @importFrom Seurat AddMetaData Cells
#' @export combineBCR
combineBCR <- function(vdj, SeuratObj,
                       keep_columns = c("v_gene", "d_gene", "j_gene",
                                        "c_gene", "full_length", "productive",
                                        "cdr3", "cdr3_nt", "reads", "umis"),
                       seurat_sample_column = NULL,
                       seurat_cell_name_column = NULL)
{
  vdj$locus <- sapply(vdj$chain, function(x){
    if( x == "IGH" ) return("H") else return("L")
  })
  vdj$locus <- factor(vdj$locus, levels = c("H", "L"))
  keep_columns <- c("sample_name", "barcode", "locus", "bcr_type", keep_columns)
  vdj <- vdj[, keep_columns]
  vdj <- split(vdj, f = vdj$locus)
  colnames(vdj[[1]])[5:ncol(vdj[[1]])] <- paste0("IGH_", colnames(vdj[[1]])[5:ncol(vdj[[1]])])
  colnames(vdj[[2]])[5:ncol(vdj[[2]])] <- paste0("IGL_", colnames(vdj[[2]])[5:ncol(vdj[[2]])])
  vdj <- merge(vdj[[1]][, which(colnames(vdj[[1]]) != 'locus')],
               vdj[[2]][, which(colnames(vdj[[2]]) != 'locus')],
               by = c("sample_name", "barcode", "bcr_type"),
               all.x = TRUE, all.y = TRUE, sort = FALSE)
  # rownames(vdj) <- vdj$barcode
  cells_in_seurat <-
  if( sum(vdj$barcode %in% Cells(SeuratObj)) == 0 ){
    # something wrong with the cell names. Use the seurat_cell_name_column
    # option and merge that way
    if( ! seurat_cell_name_column %in% colnames(SeuratObj@meta.data) )
      stop(paste0("'", seurat_cell_name_column, "' cannot be found in the columns of the meta.data slot of SeuratObj."))
    merged_metadata <- SeuratObj@meta.data
    merged_metadata$rowIDs <- rownames(merged_metadata)
    merged_metadata <- merge(merged_metadata, vdj,
                             by.x = c(seurat_sample_column, seurat_cell_name_column),
                             by.y = c("sample_name", "barcode"), all.x =TRUE,
                             all.y=  FALSE, sort = FALSE)
    rownames(merged_metadata) <- merged_metadata$rowIDs# Cells(SeuratObj)
    return( Seurat::AddMetaData(
      SeuratObj,
      merged_metadata[, colnames(vdj)[which(! colnames(vdj) %in% c("sample_name", "barcode", "rowIDs"))]])
    )
  } else {
    rownames(vdj) <- vdj$barcode
    return( Seurat::AddMetaData(SeuratObj, vdj) )
  }
}

#' Add cell metadata from Seurat object into VDJ data frame
#'
#' @description
#' `AddCellMetaToVDJ` merges annotations from Seurat object (e.g. cluster annotations) into VDJ repertoire data frame.
#'
#' @details
#' The function accepts any columns found in the metadata slot of the Seurat object. You need to indicate for each column whether the data
#'
#' @param vdj a data.frame of VDJ data after running `collapseBCR`.
#' @param SeuratObj Seurat object.
#' @param metadata_col vector of column names from SeuratObj denoting the desired metadata to be extracted and added to vdj
#' @param barcode_col column name in SeuratObj metadata which holds cell barcodes. (Default: 'row.names' i.e. use row names of SeruatObj metadata as cell barcodes)
#'
#' @return VDJ data frame with additional columns taken from the Seurat object metadata slot.
#' @export AddCellMetaToVDJ
AddCellMetaToVDJ <- function(vdj, SeuratObj, metadata_col,
                             barcode_col = "row.names")
{
  metadata <- slot(SeuratObj, name = "meta.data")
  for(i in metadata_col){
    if(! i %in% colnames(metadata))
      stop(paste0("Column '", i, "' cannot be found in meta.data slot of SeuratObj."))
  }
  if( barcode_col != 'row.names' ){
    metadata_col <- c( metadata_col, barcode_col )
    metadata <- metadata[, metadata_col]
    o <- merge(vdj, metadata, by.x = "barcode", by.y = barcode_col,
               all.x = TRUE, all.y = FALSE, sort = FALSE)
  } else {
    o <- merge(vdj, metadata, by.x = "barcode", by.y = "row.names",
               all.x = TRUE, all.y = FALSE, sort = FALSE)
  }
  return(o)
}

