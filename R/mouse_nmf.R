#' Isotype sterile/productive signatures trained using mouse single-cell B cell atlas
#'
#' The breakdown of sterile/productive IgH transcripts into isotype 'signatures' using
#' non-negative matrix factorization (NMF). This object stores the NMF weights of each
#' sterile/productive IgH gene in each NMf signature.
#'
#' @format A matrix with 14 rows and 2 columns, each row corresponds to a IgH sterile/productive gene and each column its weights in the isotype signatures.
"human_nmf"
