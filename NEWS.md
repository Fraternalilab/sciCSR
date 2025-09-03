# sciCSR 0.3.3 (Sep 3, 2025)
+ Fixed problem in `getIGHmapping` returning errors when no intronic reads could be found in the supplied BAM file. Now the `junction_reads` slot of the output of `getIGHmapping` will be NULL if this is the case.
+ Fixed problem in `mergeIgHCountsToSeurat` for discrepancies between cell barcodes in Seurat object and sterile/productive transcript count matrix.

# sciCSR 0.3.2 (Jul 11, 2024)
+ moved engine to convert between Seurat (`.rds`) and Scanpy (`.h5ad`) data objects to `sceasy`.

# sciCSR 0.3.1 (Sep 12, 2023)
+ correct problems in `annotatePairing` to give correct classification of number of H/L BCR sequences annotated to each cell barcode.

# sciCSR 0.3.0 (May 31, 2023)
* new function `getIsotype` to allow grouping of cells based on any subset of the productive/sterile count matrix (i.e. allow grouping based on sterile transcript isotype)
* new function `mergIgHCountsToSeurat` adds zeros for missing cells and merge productive/sterile count matrix into Seurat object.
* unit tests implemented using `testthat`.

# sciCSR 0.2.1 (Feb 1, 2023)
* fixed issues re distributed computing in python across MacOS and Windows.
* added functions for basic scRNA-seq data preprocessing (`normalise_dimreduce`) which includes customisable pruning of variably expressed genes, collapse VDJ (or other genes) into metagenes etc (`collapseIntoMetaGenes`).

# sciCSR 0.2.0 (Jan 30, 2023)
* removed velocyto.R dependency to avoid installation problems.
* fixed minor issues for cross-platform (MacOS / Linux) installation and use.
* fixed issues with multiprocessing in TPT calculations.
* added a `NEWS.md` file to track changes to the package.
* functional vignette site via `pkgdown`.

# sciCSR 0.1.1 (Jan 16, 2023)
* added new functions for plotting TPT isotype results.
* NMF-based calculation of CSR potential.
* functions to merge repertoire.
* plot arrow-style velocity streams a-la-`scVelo`.

# sciCSR 0.1.0 (Sep 8, 2022)
* first version of working package.

