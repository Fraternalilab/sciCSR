#' Normalisation and dimensionality reduction of scRNAseq gene counts
#'
#' @description
#' \code{normalise_dimreduce} is a wrapper function around the Seurat basic data processing workflow to generate normalised gene count data, its dimensionality reduction, and derivation of cell clusters.
#'
#' @details
#' \code{normalise_dimreduce} is a wrapper function around the major basic Seurat data processing workflow and performs the following steps:
#' \itemize{
#'   \item{Calculation of \% mitochondrial transcripts}{ (\code{PercentageFeatureSet}) and subsetting to remove those beyond the cutoff given by \code{mt.percent}.}
#'   \item{Gene count normalisation}{, using either \code{SCTransform} or \code{NormalizeData}}
#'   \item{Pruning variably expressed features}{. All genes with names matching the vector of regular expression given in the argument \code{features_exclude} will be removed from this list to avoid them influencing the downstream dimensionality reduction and clustering steps. This is particularly relevant for avoiding clusters of B cells grouped by their isotypes/VDJ expression.}
#'   \item{Principal component analysis (PCA)}{ (\code{Seurat::RunPCA} function)}
#'   \item{Batch correction using \code{Harmony}}{: covariates given in \code{harmony_vars} will be removed in the Harmony regression step. (Optional, if \code{run.harmony == TRUE})}
#'   \item{UMAP dimensionality reduction}{: \code{Seurat::RunUMAP}, retaining the top principal components, each of which explain at least \code{var_explained_lim} of the variance.}
#'   \item{k-neighbor network (kNN) construction}{ (\code{Seurat::FindNeighbors})}
#'   \item{Define cell clusters based on kNN graph}{ (\code{Seurat::FindClusters})}
#' }
#'
#' @param obj Seurat object with the gene counts unnormalised.
#' @param var_explained_lim numeric, the minimum proportion of variance explained cutoff for a principal component to be included in the dimensionality reduction and clustering steps (default: 0.015, i.e. 1.5\%)
#' @param run.harmony should the package \code{Harmony} be used on the data? (Default: FALSE)
#' @param harmony_vars vector of parameters to be included in the regression step in \code{Harmony}. Variations specific to these parameters will be removed during the \code{Harmony} run.
#' @param SCT Should the \code{SCTransform} pipeline be used? If not it will follow the standard Seurat normalisation workflow (\code{NormalizeData}, \code{FindVariableFeatures}, \code{ScaleData})
#' @param mt.pattern the regular expression used to identify mitochondrial transcripts (Default: ^MT-", i.e. all gene names beginning with "MT-")
#' @param mt.percent the cutoff for mitochondrial transcript percentage, above which cells will be removed from the Seurat project as part of quality control (Default: 10, i.e. cells with more than 10\% of counts mapped to mitochondrial transcripts will be removed from the Seurat object)
#' @param features_exclude a vector of regular expressions to select genes to be IGNORED during dimensionality reduction and clustering. By default the following features were included in this list: IgH, K, L V/D/JC genes, TRA/TRB V/C genes, AC233755.1 (which encodes a V-gene-like product), IGLL, JCHAIN)
#' @param ... Arguments to be passed to various Seurat functions (\code{SCTransform}, \code{NormalizeData}, \code{FindVariableFeatures}, \code{ScaleData}, \code{RunPCA}, \code{RunUMAP}, \code{FindNeighbors}, \code{FindClusters})
#'
#' @return A Seurat object with normalised gene count data, dimensionality reduction and clustering done
#'
#' @import Seurat
#' @importFrom harmony RunHarmony
#' @export normalise_dimreduce
normalise_dimreduce <- function(obj, var_explained_lim = 0.015,
                                run.harmony = FALSE, harmony_vars = NULL,
                                SCT = FALSE, mt.pattern = "^MT-", mt.percent = 10,
                                features_exclude = c("^IGH[MDE]", "^IGHG[1-4]", "^IGHA[1-2]",
                                                     "^IG[HKL][VDJ]", "^IGKC", "^IGLC[1-7]",
                                                     "^TR[ABGD][CV]", "^AC233755.1", "^IGLL",
                                                     "^JCHAIN"), ...)
{
  if( ! inherits(obj, "Seurat") )
    stop("'obj' should be a Seurat object.")
  if( ! is.numeric( var_explained_lim) )
    stop("'var_explained_lim' should be a numeric.")
  if( length( var_explained_lim ) != 1 )
    stop("'var_explained_lim' should be a numeric of length 1.")
  if( var_explained_lim < 0 | var_explained_lim > 1)
    stop("'var_explained_lim' should be a numericb between 0 and 1")
  if( ! mt.percent | is.numeric(mt.percent) ){
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt.pattern)
    if( ! mt.percent ) mt.percent <- 10
    subset_cells <- Cells(obj)[which(obj$percent.mt < mt.percent)]
    obj <- subset(obj, cells = subset_cells)
  }
  if(SCT){
    obj <- SCTransform(obj, ...)
    assay.use <- "SCT"
  } else {
    obj <- NormalizeData(obj, ...)
    obj <- FindVariableFeatures(obj, ...)
    all.genes <- rownames(obj)
    obj <- ScaleData(obj, features = all.genes, ...)
    assay.use <- "RNA"
  }
  varfeat <- VariableFeatures(obj)
  varfeat <- varfeat[!grepl(paste(features_exclude, collapse = "|"), varfeat)]
  obj <- RunPCA(obj, features = varfeat, ...)
  if(run.harmony){
    obj <- harmony::RunHarmony(obj, assay.use = assay.use,
                               group.by.vars = harmony_vars,
                               verbose = TRUE)
  }
  var_explained <- (obj@reductions$pca@stdev)^2 / sum((obj@reductions$pca@stdev)^2 )
  dim <- max(which(var_explained > var_explained_lim))
  cat(paste0("Top ", dim, " PCs explain variance >= ",
             var_explained_lim, ".\n"))
  if(run.harmony){
    obj <- RunUMAP(obj, reduction = "harmony", dims = 1:dim, ...)
  } else {
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:dim, ...)
  }
  obj <- FindNeighbors(obj, dims = 1:dim, ...)
  obj <- FindClusters(obj, ...)
  obj
}

#' Group gene counts into metagenes to minimise individual differences
#'
#' @description
#' \code{collapseIntoMetagenes} defines metagenes which sum over the gene counts mapped to a group of genes, in a single-cell gene expression count matrix.
#'
#' @details
#' \code{collapseIntoMetagenes} can be used to group transcript counts into metagenes, to remove the effect of e.g. individual variations which leads to preference of specific genes.
#' One example is the immunoglobulin VDJ genes whose expression is specific to each B cell and is indicative of clonotype rather than cell state. By summing over all individual VDJ genes into one metagene, this avoids those individual genes to influence the downstream dimensionality projection and clustering results.
#'
#' @param countmat sparse matrix containing single-cell gene expression data. Output from \code{Seurat::Read10X} or equivalent.
#' @param metagenes_definitions a vector containing regular expressions to match gene names in the row names of \code{countmat}. For each regular expression, matched genes will be summarised into one metagene (see Details). (Default: individual metagenes for ribosomal, HLA I-major, HLA I-minor, HLA II and Ig VDJ transcripts.)
#'
#' @return a sparse matrix containing count data where all the matched genes are collapsed into metagenes with names given by the names of each element in \code{metagenes_definitions}.
#'
#' @import Matrix
#' @export collapseIntoMetagenes
collapseIntoMetagenes <- function(countmat,
                                  metagenes_definitions = c("RIBO" = "^RP[LS]|^MRP[LS]",
                                                            "HLA-Imaj" = "^HLA-[ABC]$",
                                                            "HLA-Imin" = "^HLA-[EFG]$",
                                                            "HLA-II" = "^HLA-D",
                                                            "VDJ" = "^IG[HKL][VDJ][0-9]"))
{
  metagenes <- lapply(metagenes_definitions, function(x){
    colSums(as.matrix(
      countmat[
        rownames(countmat)[grepl(x, rownames(countmat))],
      ]
    ))
  })
  names(metagenes) <- names(metagenes_definitions)
  for(item in names(metagenes_definitions)){
    n <- item
    item <- metagenes_definitions[item]
    # Remove individual genes from the matrix
    countmat <- countmat[-(which(grepl(item, rownames(countmat)))), ]
    # Append the metagene to the matrix
    countmat <- rbind(countmat, metagenes[[n]])
  }
  # Rename the added rows
  rownames(countmat[tail(1:nrow(countmat), length(metagenes_definitions))]) <- names(metagenes_definitions)
  return( countmat )
}
