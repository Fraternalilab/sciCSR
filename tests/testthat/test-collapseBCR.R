vdj <- system.file("extdata/hong_sampled_cellranger_vdj.csv", package = "sciCSR")
vdj <- read.csv(vdj, stringsAsFactors = FALSE)
# just get a subset of cells for test
barcodes <- unique(vdj$barcode)[1:100]
vdj <- vdj[which(vdj$barcode %in% barcodes), ]

test_that(desc = "collapseBCR works",
          code = {
            collapsed <- collapseBCR(vdj, format = "10X")
            expect_equal(nrow(collapsed), 181)
            expect_equal(ncol(collapsed), 34)
            expect_true("bcr_type" %in% colnames(collapsed))
            expect_equal(unique(collapsed[which(collapsed$barcode == "ACTTACTTCAGGTAAA-1"), "bcr_type"]),
                         "multi_LC_same_class")
          })


test_that(desc = "collapseBCR can generate full table",
          code = {
            collapsed <- collapseBCR(vdj, format = "10X", full.table = TRUE)
            expect_equal(length(collapsed), 2)
            expect_equal(nrow(collapsed[[1]]), 181)
            expect_equal(nrow(collapsed[[2]]), 184)
            expect_gt(nrow(collapsed[[2]]), nrow(collapsed[[1]]))
          })

test_that(desc = "format check works",
          code = {
            expect_error(collapseBCR(vdj, format = "somethingelse"),
                         "Only '10X' is allowed in the argument 'format'")
          })

test_that(desc = "fails if c_gene column cannot be found",
          code = {
            vdj2 <- vdj
            colnames(vdj2)[which(colnames(vdj2) == "c_gene")] <- "cgene"
            expect_error(collapseBCR(vdj2),
                         "undefined columns selected")
          })
