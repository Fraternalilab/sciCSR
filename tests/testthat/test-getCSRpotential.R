hong_b <- readRDS(system.file("extdata/hong_sampled_merged.rds", package = "sciCSR"))

test_that(desc = "getCSRpotential works",
          code = {
            hong_b2 <- getCSRpotential(
              SeuratObj = hong_b,
              c_gene_anno_name = "IGH_c_gene",
              reference_based = "mouse"
            )
            expect_true("isotype" %in% colnames(slot(hong_b2, "meta.data")))
            expect_true("csr_pot" %in% colnames(slot(hong_b2, "meta.data")))
            expect_equal(round(mean(hong_b2$csr_pot), 3), 0.859)
            expect_equal(sum(hong_b2$isotype == "G2b"), 344)
          })

test_that(desc = "error if wrong species is stated in reference_based",
          code = {
            expect_error(
              getCSRpotential(
                SeuratObj = hong_b,
                c_gene_anno_name = "IGH_c_gene",
                reference_based = "human"
              ), "subscript out of bounds")
          })

test_that(desc = "input checks work",
          code = {
            expect_error(
              getCSRpotential(
                SeuratObj = hong_b, ighc_count_assay_name = "IGH",
                c_gene_anno_name = "IGH_c_gene",
                reference_based = "mouse"
              ), "assay .*cannot be found in SeuatObj"
            )
            expect_error(
              getCSRpotential(
                SeuratObj = hong_b,
                c_gene_anno_name = "IGH_cgene",
                reference_based = "mouse"
              ), "column name given in 'c_gene_anno_name' is not found"
            )
            expect_error(
              getCSRpotential(
                SeuratObj = hong_b,
                c_gene_anno_name = "IGH_c_gene",
                reference_based = "mouse", mode = "high"
              ), "'mode' must be either 'furthest' or 'highest'"
            )
            expect_error(
              getCSRpotential(
                SeuratObj = hong_b,
                c_gene_anno_name = "IGH_c_gene",
                reference_based = "mouse2",
              ), "'reference_based' must be either 'human' or 'mouse'"
            )
            expect_error(
              getCSRpotential(
                SeuratObj = hong_b,
                c_gene_anno_name = "IGH_c_gene",
                reference_based = "mouse", knn_graph = NULL
              ), "'knn_graph' must either be TRUE or FALSE"
            )
            expect_error(
              getCSRpotential(
                SeuratObj = hong_b,
                c_gene_anno_name = "IGH_c_gene",
                reference_based = "mouse", ighc_slot = "count2"
              ), "named slot.*cannot be found"
            )
          })
