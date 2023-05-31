hong_b <- readRDS(system.file("extdata/hong_sampled_merged.rds", package = "sciCSR"))

test_that(desc = "getSHM works",
          code = {
            hong_b2 <- getSHM(hong_b, v_identity_anno_name = "IGH_v_identity")
            expect_true("shm" %in% colnames(slot(hong_b2, "meta.data")))
            expect_equal(round(mean(hong_b2$shm), 2), 0.14)
          })

test_that(desc = "error if v_identity_anno_name is not found",
          code = {
            expect_error(getSHM(hong_b, v_identity_anno_name = "IGH_v_identity2"),
                         "Cannot find.*in this Seurat object")
          })
