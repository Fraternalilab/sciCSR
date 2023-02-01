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

