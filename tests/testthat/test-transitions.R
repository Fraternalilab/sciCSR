hong_b <- readRDS(system.file("extdata/hong_sampled_merged.rds", package = "sciCSR"))
hong_b <- getCSRpotential(
  SeuratObj = hong_b,
  c_gene_anno_name = "IGH_c_gene",
  reference_based = "mouse"
)
hong_b <- getSHM(hong_b, v_identity_anno_name = "IGH_v_identity")

test_that(desc = "convertSeuratToH5ad works",
          code = {
            convertSeuratToH5ad(hong_b, assays = c("RNA"), "hong_sampled_bcells.h5ad")
            expect_true(file.exists("hong_sampled_bcells_assay-RNA.h5ad"))
          })

test_that(desc = "splitAnnData works",
          code = {
            splitAnnData(anndata_file = "hong_sampled_bcells_assay-RNA.h5ad",
                         split.by = "Status", levels = c("WT", "IL23-/-"))
            expect_true(file.exists("hong_sampled_bcells_assay-RNA_IL23--.h5ad"))
          })

test_that(desc = "fitTransitionModel works",
          code = {
            g <- fitTransitionModel(
              anndata_file = "hong_sampled_bcells_assay-RNA_IL23--.h5ad",
              mode = "pseudotime", pseudotime_key = "csr_pot"
            )
            expect_equal(length(g), 3)
            expect_equal(names(g), c("cellrank_obj", "transition_matrix", "CellID"), ignore_attr = TRUE)
            expect_equal(class(g$transition_matrix), c("matrix", "array"), ignore_attr = TRUE)
            expect_equal(round(g$transition_matrix[292, 749], 3), 0.165)
            expect_equal(round(g$transition_matrix[292, 526], 3), 0.008)
            expect_equal(round(max(g$transition_matrix), 3), 0.416)
          })

test_that(desc = "fitTransitionModel works",
          code = {
            g <- fitTransitionModel(
              anndata_file = "hong_sampled_bcells_assay-RNA_IL23--.h5ad",
              mode = "pseudotime", pseudotime_key = "csr_pot"
            )
            expect_equal(length(g), 3)
            expect_equal(names(g), c("cellrank_obj", "transition_matrix", "CellID"), ignore_attr = TRUE)
            expect_equal(class(g$transition_matrix), c("matrix", "array"), ignore_attr = TRUE)
            expect_equal(round(g$transition_matrix[292, 749], 3), 0.165)
            expect_equal(round(g$transition_matrix[292, 526], 3), 0.008)
            expect_equal(round(max(g$transition_matrix), 3), 0.416)
          })

test_that(desc = "fitTPT works",
          code = {
            g <- fitTransitionModel(
              anndata_file = "hong_sampled_bcells_assay-RNA_IL23--.h5ad",
              mode = "pseudotime", pseudotime_key = "csr_pot"
            )
            tpt <- fitTPT(
              anndata_file = "hong_sampled_bcells_assay-RNA_IL23--.h5ad",
              CellrankObj = g, group.cells.by = "isotype",
              source_state = 'M', target_state = 'A', random_n = 3
            )
            expect_equal(length(tpt), 11)
            expect_equal(round(tpt$total_gross_flux, 3), 0.003)
            expect_named(tpt$stationary_distribution, expected = c("A", "G1", "G2b", "G2c", "G3", "M"))
            expect_equal(unname(round(tpt$stationary_distribution["G1"], 3)), 0.386)
            expect_equal(round(tpt$gross_flux["M", "G1"], 3), 9.724)
            expect_equal(round(tpt$significance["M", "G1"], 3), 0.064)
          })

