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
  if( class( bam_entries ) == 'try-error' ){
    # try edit the seqlevel style (i.e. '1' --> 'chr1')
    GenomeInfoDb::seqlevelsStyle( gRange ) <- 'UCSC'
    params <- Rsamtools::ScanBamParam(which = gRange,
                                      what = c("qname", "rname", "pos", "cigar", "flag"),#, "seq", "qual"),
                                      tag = tags,
                                      reverseComplement = FALSE,
                                      flag = flags)
    bam_entries <- try( GenomicAlignments::readGAlignments(file = bam, index = paste0(bam, '.bai'), 
                                                           param = params), silent = TRUE )
    if( class( bam_entries ) == 'try-error'){
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


getIGHreadType <- function(tb)
{
  read_info <- tb[, c("CB", "UB")]
  tb <- as.matrix( tb[, -which( colnames( tb ) %in% c("CB", "UB") )] )
  columns <- colnames( tb )
  c_col <- which(grepl("C$", columns))
  i_col <- which(grepl("I$", columns))
  j_col <- which(columns == "J")
  out <- apply(tb, MARGIN = 1, function(x){
    o <- c("isotype" = NA, "transcript_type" = NA)
    # first check whether it is a M/D/G/A/E only read
    # if not, indeterminate. set NA_NA
    c_reads <- x[c_col]; names(c_reads) <- columns[c_col]
    i_reads <- x[i_col]; names(i_reads) <- columns[i_col]
    all_c <- stringr::str_extract_all( names(c_reads[ c_reads > 0]), "IGH[MDEGA]" )
    all_i <- stringr::str_extract_all( names(i_reads[ i_reads > 0]), "IGH[MDEGA]" )
    all_c <- unique(unlist(all_c))
    all_i <- unique(unlist(all_i))
    if( length(all_c) > 1 | length(all_i) > 1 ){
      o <- c("isotype" = NA, "transcript_type" = NA)
      return( paste(o, collapse = "_") )
    }
    # check the intronic reads first
    if( all(i_reads == 0) ){
      # if all 0, it is not a GLT - either JC or uninformative
      # check majority-voted C gene on the CDS reads and report
      # based on whether J reads are found
      c_isotypes <- gsub("_C$", "", columns[ c_col[which(c_reads == max(c_reads))] ], "")
      if( length(c_isotypes) == 1 ){
        if( x[j_col] > 0 ) o <- c("isotype" = c_isotypes, "transcript_type" = "JC")
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
      # (1) If no CDS reads at all, either set GLT or JC depending on J reads
      #     (if unique intronic isotype) or set indeterminate (since failed to
      #     resolve isotype using BOTH I and C reads)
      if( all( c_reads == 0 ) ){
        if( length(i_isotypes) == 1 ){
          if( x[j_col] > 0 ) o <- c("isotype" = i_isotypes, "transcript_type" = "JC")
          else o <- c("isotype" = i_isotypes, "transcript_type" = "IC")
          return( paste(o, collapse = "_") )
        } else {
          o <- c("isotype" = NA, "transcript_type" = NA)
          return( paste(o, collapse = "_") )
        }
      }
      # check the CDS reads; get the majority-voted C gene
      c_isotypes <- gsub("_C$", "", columns[ c_col[which(c_reads == max(c_reads))] ], "")
      if( length(i_isotypes) == 1 ){
        # (2) check re JC or GLT depending on whether CDS reads of the 
        #     *same isotype* are found
        # (3) If CDS reads are found but isotype does not agree with the 
        #     majority-voted intronic isotype --> indeterminate. set NA_NA
        # (4) If CDS reads are found and isotype_i and _c agrees, determine GLT or JC
        if( length(c_isotypes) >= 1 && i_isotypes %in% c_isotypes ){
          if( x[j_col] > 0 ) o <- c("isotype" = i_isotypes, "transcript_type" = "JC")
          else o <- c("isotype" = i_isotypes, "transcript_type" = "IC")
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
          if( x[j_col] > 0 ) o <- c("isotype" = c_isotypes, "transcript_type" = "JC")
          else o <- c("isotype" = c_isotypes, "transcript_type" = "IC")
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

summariseIGHreads <- function(tb, IGHC_granges)
{
  tb$anno <- factor(tb$anno, 
                    levels = c(unlist(lapply(names(IGHC_granges), 
                                           function(x) paste(x, c("IC", "JC", "C"), sep = "_")))))
  reshape2::acast(plyr::ddply(tb, c("CB", "anno"), nrow), CB ~ anno, 
                  fill = 0, value.var = "V1", drop = FALSE)
#  list("GLT" = reshape2::acast(plyr::ddply(tb[grepl("GLT$", tb$anno), ], 
#                                           c("CB", "anno"), nrow), CB ~ anno, 
#                               fill = 0, value.var = "V1"),
#       "switched" = reshape2::acast(plyr::ddply(tb[grepl("switched$", tb$anno), ], 
#                                           c("CB", "anno"), nrow), CB ~ anno, 
#                               fill = 0, value.var = "V1"),
#       "uninformative" = reshape2::acast(plyr::ddply(tb[grepl("uninformative$", tb$anno), ], 
#                                           c("CB", "anno"), nrow), CB ~ anno, 
#                               fill = 0, value.var = "V1")
#       )
}

getIGHmapping <- function(bam, IGHC_granges, IGHVDJ_granges, 
                          cellBarcodeTag = "CB", umiTag = "UB",
                          paired = FALSE, flank_size = 5000)
{
  message("Fetching reads mapped to VDJ genes ...\n")
  J <- scanBam(bam, range(IGHVDJ_granges), cellBarcodeTag = cellBarcodeTag,
               umiTag = umiTag, paired = paired)
  if( nrow(J) > 0 ) J$type <- "J"
  message("Fetching reads mapped to C gene coding regions ...\n")
  C <- do.call("rbind", lapply(names(IGHC_granges), function(isotype){
    out <- scanBam(bam, IGHC[isotype], cellBarcodeTag = cellBarcodeTag,
                   umiTag = umiTag, paired = paired)
    #out <- unique(out[, c("CB", "UB")])
    if( nrow(out) > 0 ) out$type <- paste0(isotype, "_C")
    out
  }))
  message("Fetching reads mapped to C gene 5' regions ...\n")
  intronic <- do.call("rbind", lapply(names(IGHC_granges), function(isotype){
    out <- scanBam(bam, GenomicRanges::flank(IGHC[isotype], flank_size, start = FALSE), 
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
            "J")
  if( nrow(o) == 0 ) return( data.frame() )
  o <- reshape2::dcast(o, as.formula( paste0( cellBarcodeTag, " + ",
                                              umiTag, " ~ type")), 
                       value.var = "V1", fill = 0)
  for( col in cols[ which( !cols %in% colnames( o ) )]){
    o[, col] <- 0
  }
  return( list( "read_count" = o[, cols], "junction_reads" = junction_reads ) )
}

getJunctionReads <- function(read_info)
{
  # reads (& molecules) which capture different C genes (either CDS and/or Intronic)
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
