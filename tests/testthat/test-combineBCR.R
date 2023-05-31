hong_sampled <- system.file("extdata/hong_sampled.rds", package = "sciCSR")
hong_sampled <- readRDS(hong_sampled)

vdj <- system.file("extdata/hong_sampled_cellranger_vdj.csv", package = "sciCSR")
vdj <- read.csv(vdj, stringsAsFactors = FALSE)
collapsed <- collapseBCR(vdj)

test_that(desc = "error if sample_name column not in the collapsed vdj data.frame",
          code = {
            expect_error(
              combineBCR(
                collapsed, hong_sampled,
                # list the columns from vdj you wish to add to the Seurat object down here
                keep_columns = c("v_gene", "d_gene", "j_gene", "c_gene",
                                 "full_length",
                                 "productive", "cdr3", "cdr3_nt",
                                 "reads", "umis")
              ), "undefined columns selected"
            )
          })

test_that(desc = "error if column names listed in keep_columns not in the collapsed vdj data.frame",
          code = {
            colnames(collapsed)[which(colnames(collapsed) == "donor")] <- "sample_name"
            expect_error(
              combineBCR(
                collapsed, hong_sampled,
                # list the columns from vdj you wish to add to the Seurat object down here
                keep_columns = c("v_gene", "d_gene", "j_gene", "c_gene",
                                 "v_identity", "full_length",
                                 "productive", "cdr3", "cdr3_nt",
                                 "reads", "umis")
              ), "undefined columns selected"
            )
          })

test_that(desc = "output checks",
          code = {
            colnames(collapsed)[which(colnames(collapsed) == "donor")] <- "sample_name"
            hong_sampled2 <- combineBCR(
                collapsed, hong_sampled,
                # list the columns from vdj you wish to add to the Seurat object down here
                keep_columns = c("v_gene", "d_gene", "j_gene", "c_gene",
                                 "full_length",
                                 "productive", "cdr3", "cdr3_nt",
                                 "reads", "umis")
            )
            expect_equal(ncol(slot(hong_sampled2, "meta.data")), 31)
            expect_true("IGH_productive" %in% colnames(slot(hong_sampled2, "meta.data")))
            expect_true("IGL_c_gene" %in% colnames(slot(hong_sampled2, "meta.data")))
            expect_true("bcr_type" %in% colnames(slot(hong_sampled2, "meta.data")))
          })
