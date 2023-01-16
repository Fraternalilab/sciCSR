#' Scan reads mapped to a given genomic range
#'
#' @description
#' `scanBam` looks for reads in a BAM file which are mapped to a genomic
#' range given by the parameter `gRange`.
#'
#' @details
#' The function looks for BAM alignments which are mapped to the genomic
#' range given by the parameter `gRange`. All alignments which *span across*
#' `gRange` but do not have bases which reside within the range will be ignored
#' (e.g. a spliced read spans across a given range but will not have any bases which
#' map to the spliced range). The function expects [standard BAM/SAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
#' plus additional tags indicating the associated cell and molecule barcodes.
#' Set `cellBarcodeTag` and/or `umiTag` as NULL or NA if these tags are absent from
#' the given BAM file.
#'
#' @param bam filepath to the BAM file to read
#' @param gRange `GenomicRanges::GRanges` object specifying the genomic range to scan. See Examples.
#' @param cellBarcodeTag Name of tag holding information about cell barcode. The code expects and extracts this tag from each line in the BAM alignments. Set as NULL or NA if no such information is available in the BAM file. (Default: "CB")
#' @param umiTag Name of tag holding information about molecule barcode. The code expects and extracts this tag from each line in the BAM alignments. Set as NULL or NA if no such information is available in the BAM file. (Default: "UB")
#' @param paired Are the sequencing reads paired-end? (Default: FALSE)
#'
#' @return A data.frame:
#' \describe{
#'   \item{qname}{read ID in the FASTQ/BAM}
#'   \item{rname}{chromosome}
#'   \item{pos}{genomic position}
#'   \item{cigar}{CIGAR string indicating the result of aligning the sequencing read to the nominated position in the genome. See the SAM format specification document linked above in the Description for details}
#'   \item{flag}{the FLAG field in the SAM file. Indication of the nature of the alignment. See the SAM specification linked above in the Description for more details.}
#'   \item{`cellBarcodeTag` (Default: "CB")}{cell barcode}
#'   \item{`umiTag` (Default: "UB")}{Unique molecule identifier (UMI)}
#' }
#' The scanning of the BAM file by default removes secondary alignments, duplicates and reverse complement alignments.
#'
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom GenomicAlignments readGAlignments findOverlaps
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomicRanges mcols
#' @importFrom S4Vectors countQueryHits
#'
#' @export scanBam

scanBam <- function(bam, gRange, cellBarcodeTag = "CB", umiTag = "UB",
                    paired = FALSE)
{
  tags <- c(cellBarcodeTag, umiTag)
  tags <- tags[!is.null(tags)]
  tags <- tags[!is.na(tags)]
  if( length( tags ) == 0 ) tags <- character(0)
  flags = Rsamtools::scanBamFlag(isSecondaryAlignment=NA,#FALSE,
                                 isDuplicate=NA)#FALSE)
  if( paired ){
    flags <- Rsamtools::scanBamFlag(isSecondaryAlignment=NA,#FALSE,
                                    isDuplicate=NA,#FALSE,
                                    isPaired=TRUE,
                                    isProperPair=TRUE, isFirstMateRead=TRUE)
  }
  params <- Rsamtools::ScanBamParam(which = gRange,
                                    what = c("qname", "rname", "pos", "cigar", "flag"),#, "seq", "qual"),
                                    tag = tags,
                                    reverseComplement = FALSE,
                                    flag = flags)
  bam_entries <- list()
  if ( !file.exists(  paste0(bam, '.bai') ) ){
    stop("BAM index (*.bai) file is not present. Sort your bam file on the command line with e.g. samtools sort <BAM_FILE> -o <BAM_SORTED> if you haven't done so, and then generate BAM index on the command line with: samtools index <BAM_SORTED>")
  }
  bam_entries <- try( GenomicAlignments::readGAlignments(file = bam, index = paste0(bam, '.bai'),
                                                         param = params), silent = TRUE )
  if( inherits( bam_entries, 'try-error' ) ){
    # try edit the seqlevel style (i.e. '1' --> 'chr1')
    GenomeInfoDb::seqlevelsStyle( gRange ) <- 'UCSC'
    params <- Rsamtools::ScanBamParam(which = gRange,
                                      what = c("qname", "rname", "pos", "cigar", "flag"),#, "seq", "qual"),
                                      tag = tags,
                                      reverseComplement = FALSE,
                                      flag = flags)
    bam_entries <- try( GenomicAlignments::readGAlignments(file = bam, index = paste0(bam, '.bai'),
                                                           param = params), silent = TRUE )
    if( inherits( bam_entries, 'try-error') ){
      stop("Failed to fetch any alignment in the given genomic region. Are you sure you have indicated the correct BAM file and/or genomic region?")
    }
  }
  # filter away reads which *span* the given gRange but doesn't have matches in the range
  in_range <- GenomicAlignments::findOverlaps(bam_entries, gRange)
  GenomicRanges::mcols(bam_entries)$within <- S4Vectors::countQueryHits(in_range)
  out <- as.data.frame(GenomicRanges::mcols(bam_entries))
  out <- out[out$within == 1, ]
  if( !is.null( cellBarcodeTag ) && !is.na( cellBarcodeTag )) {
    out <- out[!is.na(out[, cellBarcodeTag]), ]
  }
  if( !is.null( umiTag) && !is.na( umiTag )) {
    out <- out[!is.na(out[, umiTag]), ]
  }
  unique(out)
}

#' Deduce type of IGH reads
#'
#' @description
#' `getIGHreadType` deduces, for each molecule in each cell identified via `scanBam`,
#' whether it corresponds to a productive (labelled "-P") or a sterile ("-S") IGH
#' molecule. A third category ("-C") is assigned if insufficient information is present
#' to distinguish "-P" and "-S".
#'
#' @details
#' The function loops through each cell barcode - molecule barcode combination
#' and looks for evidence of reads mapping to the following three regions: VDJ,
#' C as well as the 5' intron to C. Each molecule is assigned based on these criteria:
#' \itemize{
#'  \item{Productive (-P)} {at least one read mapping to VDJ region and at least one
#' read mapping to the C exonic region. NOTE: here we are not making the distinction
#' of whether the VDJ has already been spliced to the C exons at the RNA level. We are
#' simply asking whether the molecule can encode a Ig protein with both V and C regions.}
#'  \item{Sterile (-S)} {at least one read mapping to the 5' intronic region to a C gene,
#'  without reads mapping to the VDJ region. These represent the 'sterile'/'germline'
#'  IgH transcripts which primes class-switch recombination.}
#' }
#' A third category (-C) is assigned if insufficient information is available to classify
#' the molecule in to the '-P' or the '-S' groups.
#'
#' @param tb A data.frame, output from `scanBam`.
#'
#' @importFrom stringr str_extract_all
#' @export getIGHreadType

getIGHreadType <- function(tb)
{
  read_info <- tb[, c("CB", "UB")]
  tb <- as.matrix( tb[, -which( colnames( tb ) %in% c("CB", "UB") )] )
  columns <- colnames( tb )
  c_col <- which(grepl("C$", columns))
  i_col <- which(grepl("I$", columns))
  j_col <- which(columns == "VDJ")
  out <- apply(tb, MARGIN = 1, function(x){
    o <- c("isotype" = NA, "transcript_type" = NA)
    # first check whether it is a M/D/G/A/E only read
    # if not, indeterminate. set NA_NA
    c_reads <- x[c_col]; names(c_reads) <- columns[c_col]
    i_reads <- x[i_col]; names(i_reads) <- columns[i_col]
    all_c <- stringr::str_extract_all( names(c_reads[ c_reads > 0]), "I[Gg][Hh][MmDdEeGgAaTtYyXxWw]" )
    all_i <- stringr::str_extract_all( names(i_reads[ i_reads > 0]), "I[Gg][Hh][MmDdEeGgAaTtYyXxWw]" )
    all_c <- unique(unlist(all_c))
    all_i <- unique(unlist(all_i))
    if( length(all_c) > 1 | length(all_i) > 1 ){
      o <- c("isotype" = NA, "transcript_type" = NA)
      return( paste(o, collapse = "_") )
    }
    # check the intronic reads first
    if( all(i_reads == 0) ){
      # if all 0, it is not a Sterile - either Productive or uninformative
      # check majority-voted C gene on the CDS reads and report
      # based on whether J reads are found
      c_isotypes <- gsub("_C$", "", columns[ c_col[which(c_reads == max(c_reads))] ], "")
      if( length(c_isotypes) == 1 ){
        if( x[j_col] > 0 ) o <- c("isotype" = c_isotypes, "transcript_type" = "P")
        else o <- c("isotype" = c_isotypes, "transcript_type" = "C")
        return( paste(o, collapse = "_") )
      } else {
        # indeterminate. Set NA_NA
        o <- c("isotype" = NA, "transcript_type" = NA)
        return( paste(o, collapse = "_") )
      }
    } else {
      # check the intronic reads; get the majority-voted C gene
      i_isotypes <- gsub("_I$", "", columns[ i_col[which(i_reads == max(i_reads))] ], "")
      # majority voting on CDS reads:
      # (1) If no CDS reads at all, either set Sterile or Productive depending on J reads
      #     (if unique intronic isotype) or set indeterminate (since failed to
      #     resolve isotype using BOTH I and C reads)
      if( all( c_reads == 0 ) ){
        if( length(i_isotypes) == 1 ){
          if( x[j_col] > 0 ) o <- c("isotype" = i_isotypes, "transcript_type" = "P")
          else o <- c("isotype" = i_isotypes, "transcript_type" = "S")
          return( paste(o, collapse = "_") )
        } else {
          o <- c("isotype" = NA, "transcript_type" = NA)
          return( paste(o, collapse = "_") )
        }
      }
      # check the CDS reads; get the majority-voted C gene
      c_isotypes <- gsub("_C$", "", columns[ c_col[which(c_reads == max(c_reads))] ], "")
      if( length(i_isotypes) == 1 ){
        # (2) check re Productive or Sterile depending on whether CDS reads of the
        #     *same isotype* are found
        # (3) If CDS reads are found but isotype does not agree with the
        #     majority-voted intronic isotype --> indeterminate. set NA_NA
        # (4) If CDS reads are found and isotype_i and _c agrees, determine
        #     Sterile or Productive
        if( length(c_isotypes) >= 1 && i_isotypes %in% c_isotypes ){
          if( x[j_col] > 0 ) o <- c("isotype" = i_isotypes, "transcript_type" = "P")
          else o <- c("isotype" = i_isotypes, "transcript_type" = "S")
          return( paste(o, collapse = "_") )
        } else {
          # indeterminate. Set NA_NA
          o <- c("isotype" = NA, "transcript_type" = NA)
          return( paste(o, collapse = "_") )
        }
      } else {
        # check CDS reads to determine which of the possible
        # isotypes is the most likely
        if( length(c_isotypes) == 1 && c_isotypes %in% i_isotypes ){
          if( x[j_col] > 0 ) o <- c("isotype" = c_isotypes, "transcript_type" = "P")
          else o <- c("isotype" = c_isotypes, "transcript_type" = "S")
          return( paste(o, collapse = "_") )
        } else {
          # indeterminate. Set NA_NA
          o <- c("isotype" = NA, "transcript_type" = NA)
          return( paste(o, collapse = "_") )
        }
      }
    }
  })
  out <- data.frame(read_info, out, stringsAsFactors = FALSE)
  colnames( out ) <- c("CB", "UB", "anno")
  out <- out[ which(out[, 3] != "NA_NA"), ]
  out
}

#' cast data frame IGH counts into a matrix
#'
#' @description
#' `summariseIGHreads` converts the data.frame of sterile/productive IGH molecules into a count matrix.
#'
#' @details
#' The function does a count of molecules classified to be Sterile (S), Productive (P) or C-only (C) for each isotype,
#' over molecules per cell barcode and cast this count into a matrix of cell barcodes by IGH molecule type.
#'
#' @param tb A data.frame, output from `getIGHreadType`.
#' @param IGHC_granges `GenomicRanges` object detailing the positions of IgH C genes in the genome.
#'
#' @return An array of cell barcode by IGH gene type (i.e. S/P/C per C gene)
#'
#' @importFrom reshape2 acast
#' @importFrom plyr ddply
#'
#' @export summariseIGHreads
summariseIGHreads <- function(tb, IGHC_granges)
{
  tb$anno <- factor(tb$anno,
                    levels = c(unlist(lapply(names(IGHC_granges),
                                           function(x) paste(x, c("S", "P", "C"), sep = "_")))))
  reshape2::acast(plyr::ddply(tb, c("CB", "anno"), nrow), CB ~ anno,
                  fill = 0, value.var = "V1", drop = FALSE)
}

#' wrapper function to scan sterile/productive IGH molecules from BAM file
#'
#' @description
#' `getIGHmapping` is the wrapper function intended for users to supply a BAM file,
#' and scan for sterile and productive IGH transcripts over each C gene as defined in the parameter
#' `IGHC_granges`.
#'
#' @details
#' The function reads in two `GenomicRanges::GRanges` objects, one defining the genomic
#' coordinates of the IGHC genes and another for the VDJ genes. It scans the BAM file
#' for reads covering these regions, extracting their mapped cell barcodes and Unique
#' Molecule identifier (UMI). Sterile reads are defined as those covering the intronic region
#' upstream of the 5' end of the C gene coding region - the default is to consider the region
#' (min(`previous_C_CDS_end`, `-flank`), 0), where `previous_C_CDS_end` refers to the
#' 3' end of the coding region of the previous C gene, and `flank` is an integr, given by the
#' user, which indicates 'how far' the function should looks 5' of the coding region for
#' sterile reads.
#' If you have information on where the sterile transcripts begin, these coordinatees can be
#' passed as a `GRanges` object to the parameter `flank` (see examples).
#' Outputs data frame of cell barcodes and molecules mapped to
#' each IGHC gene, classified as sterile/productive/C-only.
#'
#' @param bam filepath to the BAM file to read
#' @param IGHC_granges `GenomicRanges::GRanges` object specifying the genomic range of the IGH C gene **coding regions**.
#' @param IGHVDJ_granges `GenomicRanges::GRanges` object specifying the genomic range of the IGH VDJ genes.
#' @param cellBarcodeTag Name of tag holding information about cell barcode. The code expects and extracts this tag from each line in the BAM alignments. Set as NULL or NA if no such information is available in the BAM file. (Default: "CB")
#' @param umiTag Name of tag holding information about molecule barcode. The code expects and extracts this tag from each line in the BAM alignments. Set as NULL or NA if no such information is available in the BAM file. (Default: "UB")
#' @param paired Are the sequencing reads paired-end? (Default: FALSE)
#' @param flank either (1) an integer (indicating 5' distance from the CH exons) or (2) a `GRanges` object (indicating exact genomic positions) for defining sterile IgH transcripts. (See Examples)
#'
#' @return A list with two items:
#' \describe{
#'   \item{read_count}{a data.frame in wide format indicating for each cell barcodes and UMI combination, the number of **reads** (Note: NOT UMI!) covering VDJ, and the Coding region (C) or 5' intronic region (I) of each IGH C gene.}
#'   \item{junction_reads}{a data.frame of spliced reads and their mapped cell barcodes & UMIs. Either genuine spliced productive IgH transcripts, or strange molecules potentially worth detailed inspection.}
#' }
#'
#' @importFrom GenomicRanges sort start end flank GRanges seqnames
#' @importFrom IRanges IRanges
#' @importFrom plyr ddply
#' @importFrom reshape2 dcast
#' @importFrom stats as.formula
#'
#' @export getIGHmapping
getIGHmapping <- function(bam, IGHC_granges, IGHVDJ_granges,
                          cellBarcodeTag = "CB", umiTag = "UB",
                          paired = FALSE, flank = 5000)
{
  # first sort the coordinates of the two GRanges objs
  IGHC_granges <- GenomicRanges::sort(IGHC_granges, decreasing = TRUE)
  IGHVDJ_granges <- GenomicRanges::sort(IGHVDJ_granges, decreasing = TRUE)
  message("Fetching reads mapped to VDJ genes ...")
  J <- scanBam(bam, range(IGHVDJ_granges), cellBarcodeTag = cellBarcodeTag,
               umiTag = umiTag, paired = paired)
  if( nrow(J) > 0 ) J$type <- "VDJ"
  message("Fetching reads mapped to C gene coding regions ...")
  C <- do.call("rbind", lapply(names(IGHC_granges), function(isotype){
    out <- scanBam(bam, IGHC_granges[isotype], cellBarcodeTag = cellBarcodeTag,
                   umiTag = umiTag, paired = paired)
    #out <- unique(out[, c("CB", "UB")])
    if( nrow(out) > 0 ) out$type <- paste0(isotype, "_C")
    out
  }))
  message("Fetching reads mapped to C gene 5' regions ...")
  # pad the ranges with `flank` at the 5' end
  # in case if the padding crosses into the previous C gene, overwrite the start site
  # so it starts straight after
  if( is.numeric(flank) ){
    IGHC_flank <- c()
    # for the first gene, consider the 3' end of the last VDJ gene
    flanked <- GenomicRanges::flank(IGHC_granges[1], flank,
                                    start = FALSE) # remember the genes are on the minus strand and hence reversed in chr coords
    crossed <- (GenomicRanges::end(flanked) > min(GenomicRanges::start(IGHVDJ_granges)))
    if(crossed){
      IGHC_flank <- c(IGHC_flank, GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(IGHC_granges[1]),
                                                         ranges = IRanges::IRanges(GenomicRanges::start(IGHC_granges[1]),
                                                                                   min(GenomicRanges::start(IGHVDJ_granges)) -1)))
    } else {
      IGHC_flank <- c(IGHC_flank, GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(IGHC_granges[1]),
                                                         ranges = IRanges::IRanges(GenomicRanges::start(IGHC_granges[1]),
                                                                                   GenomicRanges::end(IGHC_granges[1]) +
                                                                                     flank)))
    }
    # for the subsequent C genes, consider the previous (ie 5') C gene
    for(i in 2:length(IGHC_granges)){
      flanked <- GenomicRanges::flank(IGHC_granges[i], flank,
                                      start = FALSE) # remember the genes are on the minus strand and hence reversed in chr coords
      crossed <- (GenomicRanges::end(flanked) > GenomicRanges::start(IGHC_granges[i - 1]))
      if( crossed ){
        IGHC_flank <- c(IGHC_flank, GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(IGHC_granges[i]),
                                                           ranges = IRanges::IRanges(GenomicRanges::start(IGHC_granges[i]),
                                                                                     GenomicRanges::start(IGHC_granges[i - 1]) -1)))
      } else IGHC_flank <- c(IGHC_flank, GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(IGHC_granges[i]),
                                                                ranges = IRanges::IRanges(GenomicRanges::start(IGHC_granges[i]),
                                                                                          GenomicRanges::end(IGHC_granges[i]) +
                                                                                            flank)))
    }
    IGHC_flank <- do.call("c", IGHC_flank)
    names(IGHC_flank) <- names(IGHC_granges)
  } else if( inherits(flank, "GRanges") ) {
    if( ! all(names(flank) == names(IGHC_granges)) ){
      stop("'flank' and 'IGHC_granges' appear to have different order of C genes. Please fix them to have same order.")
    }
    IGHC_flank <- flank
  } else {
    stop("'flank' must either be an integer (indicating 5' distance from the CH exons) or a GRanges object (indicating exact genomic positions) for defining sterile IgH transcripts.")
  }


  intronic <- do.call("rbind", lapply(names(IGHC_granges), function(isotype){
    out <- scanBam(bam, IGHC_flank[isotype],
                   cellBarcodeTag = cellBarcodeTag,
                   umiTag = umiTag, paired = paired)
    #out <- unique(out[, c("CB", "UB")])
    if( nrow(out) > 0 ) out$type <- paste0(isotype, "_I")
    out
  }))
  if( is.null( cellBarcodeTag ) || is.na( cellBarcodeTag ) ){
    cellBarcodeTag <- "CB"
    if( nrow(J) > 0 ) J[, cellBarcodeTag] <- bam # name it with the bamfile
    if( nrow(C) > 0 ) C[, cellBarcodeTag] <- bam
    if( nrow(intronic) > 0 ) intronic[, cellBarcodeTag] <- bam
  }
  if( is.null( umiTag ) || is.na( umiTag ) ){
    umiTag <- "UB"
    if( nrow(J) > 0 ) J[, umiTag] <- "None"
    if( nrow(C) > 0 ) C[, umiTag] <- "None"
    if( nrow(intronic) > 0 ) intronic[, umiTag] <- "None"
  }
  junction_reads <- getJunctionReads( rbind(J, C, intronic) )
  J <- plyr::ddply(J, c(cellBarcodeTag, umiTag, "type"), function(x) length(unique(x[, "qname"])))
  C <- plyr::ddply(C, c(cellBarcodeTag, umiTag, "type"), function(x) length(unique(x[, "qname"])))
  intronic <- plyr::ddply(intronic, c(cellBarcodeTag, umiTag, "type"),
                          function(x) length(unique(x[, "qname"])))
  o <- rbind(J, C, intronic)
  cols <- c(cellBarcodeTag, umiTag,
            unlist(lapply(names(IGHC_granges),
                          function(x) paste(x, c("C", "I"), sep = "_"))),
            "VDJ")
  if( nrow(o) == 0 ) return( data.frame() )
  o <- reshape2::dcast(o, as.formula( paste0( cellBarcodeTag, " + ",
                                              umiTag, " ~ type")),
                       value.var = "V1", fill = 0)
  for( col in cols[ which( !cols %in% colnames( o ) )]){
    o[, col] <- 0
  }
  return( list( "read_count" = o[, cols], "junction_reads" = junction_reads ) )
}

#' extract reads covering splice junctions
#'
#' @description
#' `getJunctionReads` extract reads from `scanBam` that cover splice junctions.
#'
#' @details
#' `getJunctionReads` extract reads from `scanBam` that cover splice junctions.
#' This considers only those which cover BOTH VDJ and C exons/introns, and
#' strange cases which cover multiple different IGHC genes.
#'
#' @param read_info output from `scanBam` containing details of the aligned reads for VDJ, C and intronic regions.
#'
#' @return A data.frame containing details of reads which cover splice junctions:
#' \describe{
#'   \item{qname}{read ID in the FASTQ/BAM}
#'   \item{CB)}{cell barcode}
#'   \item{UB)}{Unique molecule identifier (UMI)}
#'   \item{type}{a list of mapped entities (VDJ/IGHC genes/intronic rgeions) found on the read.}
#' }
#'
#' @importFrom plyr ddply
getJunctionReads <- function(read_info)
{
  # extract the spliced reads (ie those which capture the splice junctions)
  # return data.frame of CB, UB, qname and list of junctions captured

  # spliced reads (i.e. those with 'N' in cigar)
  read_info <- read_info[ which(grepl("N", read_info$cigar)), ]
  # eliminate splicing within a C gene CDS by selecting only those appearing in
  # the scans under different C genes
  multiple_genes <- read_info[, c("qname", "CB", "UB", "flag", "type")]
  # multiple_genes$type <- gsub("_C$|_I$", "", multiple_genes$type)
  multiple_genes <- unique( multiple_genes )
  dups <- unique( multiple_genes[duplicated(multiple_genes[, 1:4]), 1:4] )
  dups$duplicate <- TRUE
  reads <- merge( multiple_genes, dups, all.x = TRUE, sort = FALSE)
  reads <- reads[ which( !is.na( reads$duplicate ) ), ]
  reads <- plyr::ddply( reads, 1:4, function(x) paste(sort(x[, 5]), collapse = "-"))
  reads <- plyr::ddply( reads, 1:3, function(x) paste(unique(x[, 5]), collapse = ";"))
  # disregard reads concerning CDS of
  reads
}
