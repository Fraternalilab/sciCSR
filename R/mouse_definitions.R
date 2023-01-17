#' Genomic coordinates of mouse heavy-chain immunoglobulin V, D, J, C genes
#'
#' Downloaded from Ensembl BioMart query of the mm10 reference genome. The positions of sterile C transcripts were separately annotated. This is used for enumerating the productive and sterile IgH transcripts from given BAM files.
#'
#' @format A list of three GenomicRanges objects:
#' \describe{
#'   \item{VDJ}{GenomicRanges object containing genomic coordinates of individual V, D and J genes}
#'   \item{C}{GenomicRanges object containing genomic coordinates of the coding segment of individual C genes}
#'   \item{sterile}{GenomicRanges object containing genomic coordinates of the sterile transcripts of individual C genes}
#' }
"mouse_definitions"
