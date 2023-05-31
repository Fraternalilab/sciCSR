data("mouse_definitions")

hong_b <- system.file("extdata/hong_sampled.rds", package = "sciCSR")
hong_b <- readRDS(hong_b)

bamfiles <- c(
  system.file("extdata/Hong_S1_sampled_Igh.bam", package = "sciCSR"),
  system.file("extdata/Hong_S2_sampled_Igh.bam", package = "sciCSR"),
  system.file("extdata/Hong_S3_sampled_Igh.bam", package = "sciCSR"),
  system.file("extdata/Hong_S4_sampled_Igh.bam", package = "sciCSR"),
  system.file("extdata/Hong_S5_sampled_Igh.bam", package = "sciCSR")
)

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

test_that(desc = "check things work",
          code = {
            hong_IGH2 <- repairBarcode(
              hong_IGH, hong_b,
              sample_names = c("p19kd_1", "p19kd_2", "p19kd_3", "WT_4", "WT_5"),
              seurat_sample_column = "Sample"
            )
            expect_identical(unname(rownames(hong_IGH2[[2]])[34]), "ACGGGCTTCATCTGTT-2")
            expect_identical(rownames(hong_IGH[[5]])[346], "GAATGAAGTAAGCACG-1")
            expect_identical(unname(rownames(hong_IGH2[[5]])[346]), "GAATGAAGTAAGCACG-5")
            expect_match(rownames(hong_IGH2[[4]])[148], "AGAGCGATCTTCATGT")
          })

test_that(desc = "error when wrong seurat_sample_column supplied",
          code = {
            expect_error(
              repairBarcode(
                hong_IGH, hong_b,
                sample_names = c("p19kd_1", "p19kd_2", "p19kd_3", "WT_4", "WT_5"),
                seurat_sample_column = "Status"
              ),
              "incorrect number of dimensions"
            )
          })

