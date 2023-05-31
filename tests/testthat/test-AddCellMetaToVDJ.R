hong_sampled <- system.file("extdata/hong_sampled.rds", package = "sciCSR")
hong_sampled <- readRDS(hong_sampled)

vdj <- system.file("extdata/hong_sampled_cellranger_vdj.csv", package = "sciCSR")
vdj <- read.csv(vdj, stringsAsFactors = FALSE)
collapsed <- collapseBCR(vdj)

test_that(desc = "AddCellMetaToVDJ works",
          code = {
            vdj2 <- AddCellMetaToVDJ(vdj, hong_sampled,
                                     metadata_col = c("Sample", "Status", "percent.mt"),
                                     barcode_col = "row.names")
            expect_equal(ncol(vdj2), ncol(vdj) + 3)
            expect_equal(nrow(vdj), nrow(vdj2))
            expect_equal(vdj2[123, "Status"], "IL23-/-")
            expect_equal(round(vdj2[4567, "percent.mt"], 2), 2.76)
          })

test_that(desc = "checks for columns existing in Seurat object work",
          code = {
            expect_error(
              AddCellMetaToVDJ(vdj, hong_sampled,
                               metadata_col = c("Sample", "Status", "percentmt"),
                               barcode_col = "row.names"),
              "Column 'percentmt' cannot be found in meta.data slot of SeuratObj"
            )
          })

