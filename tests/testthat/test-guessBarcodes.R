test_that(desc = "guessBarcodes works on 10X-styled barcodes",
          code = {
            out <- guessBarcodes("ATCGATGGATAGCCTA-13")
            expect_true(is.na(out[1]))
            expect_equal(out[2], "ATCGATGGATAGCCTA")
            expect_equal(out[3], "-13")
          })

test_that(desc = "guessBarcodes reject short barcodes",
          code = {
            expect_error(guessBarcodes("ATCGA-13"), "the given cell_name doesn't appear to contain nucleotide barcode strings")
          })

test_that(desc = "error if min_barcode_length is not integer",
          code = {
            expect_error(guessBarcodes("ATCGATGGATAGCCTA-13", 2.5), "min_barcode_length should be an integer")
          })

test_that(desc = "error if non nucleotide-style cell barcode is supplied",
          code = {
            expect_error(guessBarcodes("Cell_1234567"), "the given cell_name doesn't appear to contain nucleotide barcode strings")
          })

