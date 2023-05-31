data("mouse_definitions")

bamfiles <- c(
  system.file("extdata/Hong_S1_sampled_Igh.bam", package = "sciCSR"),
  system.file("extdata/Hong_S2_sampled_Igh.bam", package = "sciCSR"),
  system.file("extdata/Hong_S3_sampled_Igh.bam", package = "sciCSR"),
  system.file("extdata/Hong_S4_sampled_Igh.bam", package = "sciCSR"),
  system.file("extdata/Hong_S5_sampled_Igh.bam", package = "sciCSR")
)

test_that(desc = "getIGHmapping works",
          code = {
            out <- getIGHmapping(bamfiles[1], mouse_definitions)
            expect_equal(length(out), 2)
            expect_equal(names(out), c("read_count", "junction_reads"), ignore_attr = TRUE)
            expect_equal(out$read_count[2, "Ighg2c_C"], 2)
            expect_equal(out$read_count[5, "VDJ"], 2)
            expect_equal(out$read_count[453, "Ighg2b_I"], 1)
            expect_equal(out$read_count[885, "VDJ"], 2)
            expect_equal(out$read_count[885, "Ighg2b_C"], 1)
          })

test_that(desc = "getIGHreadType works",
          code = {
            out <- getIGHmapping(bamfiles[1], mouse_definitions)
            out2 <- getIGHreadType(out$read_count)
            expect_equal(ncol(out2), 3)
            expect_equal(nrow(out2), 7009)
            expect_equal(out2[which(out2$CB == "ATCATCTTCGTTTAGG-1" &
                                      out2$UB == "CATTTTCCAG"), "anno"], "Ighg2c_P")
            expect_equal(out2[which(out2$CB == "CAGTCCTTCTCCGGTT-1" &
                                      out2$UB == "CTCTACAGTG"), "anno"], "Igha_S")
          })

test_that(desc = "summariseIGHreads works",
          code = {
            out <- getIGHmapping(bamfiles[1], mouse_definitions)
            out2 <- getIGHreadType(out$read_count)
            out3 <- summariseIGHreads(out2, mouse_definitions)
            expect_equal(dim(out3), c(519, 21), ignore_attr = TRUE)
            expect_equal(out3["AAACGGGGTAAGGGCT-1", "Ighm_C"], 11)
            expect_equal(out3["AAAGATGCAAACAACA-1", "Ighg3_S"], 1)
            expect_equal(out3["AACTCAGTCATTATCC-1", "Ighg2b_S"], 19)
            expect_equal(out3["AACTCAGTCATTATCC-1", "Ighm_P"], 2)
          })

# loop through these BAM files
hong_IGH <- lapply(bamfiles, function(bamfile){
  cat(paste0(bamfile, " ...\n"))
  # These three functions count productive and sterile
  # reads from the BAM file and give a count matrix
  out <- getIGHmapping(bamfile, mouse_definitions)
  out2 <- getIGHreadType(out$read_count)
  out3 <- summariseIGHreads(out2, mouse_definitions)
  # 'out3' is a count matrix of productive/sterile transcript
  # of each isotype per cell
  out3
})
hong_sampled <- readRDS(system.file("extdata/hong_sampled.rds", package = "sciCSR"))
hong_IGH2 <- repairBarcode(
  hong_IGH, hong_sampled,
  sample_names = c("p19kd_1", "p19kd_2", "p19kd_3", "WT_4", "WT_5"),
  seurat_sample_column = "Sample"
)

# combine these individual count matrices
hong_IGH2 <- do.call("rbind", hong_IGH2)

# remove cells which are not in the Seurat object but happen to have
# observed productive/sterile transcripts
hong_IGH2 <- hong_IGH2[which(rownames(hong_IGH2) %in% Seurat::Cells(hong_sampled)), ]

test_that(desc = "mergeIgHCountsToSeurat works",
          code = {
            hong_combined2 <- mergeIgHCountsToSeurat(
              igh_counts = hong_IGH2, SeuratObj = hong_sampled,
              assay_name = "IGHC"
            )
            expect_equal(Seurat::Assays(hong_combined2), c("RNA", "IGHC"), ignore_attr = TRUE)
            expect_equal(dim(slot(hong_combined2, "assays")[["IGHC"]]), c(21, 3000), ignore_attr = TRUE)
          })
