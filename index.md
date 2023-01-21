# sciCSR (Single-cell Inference of Class-Switch Recombination)

## Introduction

`sciCSR` is a R package designed to analyse class-switch recombination (CSR) in B cell single-cell RNA sequencing (scRNA-seq) data. CSR is a major process in B cell maturation whereby B cell changes the constant region ("isotype") of their B cell receptor (BCR) in order to adapt their function to different tissue and biological contexts. Along with somatic hypermutation (SHM) where mutations are accumulated in the variable regions of BCR, CSR is a hallmark of the transitions of B cells from a naive state to acquire memory against antigens. In sciCSR we provides routines to extract information on CSR and SHM from scRNA-seq data and, if available, scBCR-seq (i.e. single-cell BCR repertoire profiling) data.

In theory CSR and SHM information can be used as alternatives to the popular [RNA velocity]() method to infer transitions between B cells. sciCSR takes forward the extracted CSR/SHM information and inputs them to the [CellRank]() method to infer transitions between B cells; this can be used directly to analyse patterns of B cell maturation in the data, by considering evidence from CSR and SHM, both of which are native to B cell biology, addressing the [limitations]() of RNA velocity in analysing scRNA-seq datasets of mature cell types such as B cells found in circulation or in secondary lymphoid organs.

## Installation

```
# install dependency 'SeuratDisk' (requires 'remotes')
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

# install dependencies and sciCSR itself (requires 'devtools')
if (!require('devtools')){
  install.packages('devtools')
}
devtools::install_github("velocyto-team/velocyto.R")
devtools::install_github('Fraternalilab/sciCSR')
```

sciCSR requires a few python packages to be installed ([scanpy](), [scvelo](), cellrank etc.) as its functionalities depend on these python packages. We recommend setting up a [conda environment]() with these dependencies so that sciCSR exclusively calls packages from this environment; this will avoid problems with messing up your local installation of python (if any) if you have previously set up these packages for other uses. **Before the first time you run any examples/analyses with sciCSR, run the following lines in R** after finishing the sciCSR installation:

```{r}
library(sciCSR)
prepare_sciCSR()
```

This will set up an environment named 'scicsr' with all the desired dependencies ready for your analysis. The code should detect any existing Anaconda/Miniconda installation and uses that to set up the environment; if not, it will install Miniconda at the R default location and set up the conda environment.

## Examples

Please consult the following vignettes:

* [**Analysing CSR**](): This vignette analyses a subset of data from [Hong et al. J Immunol 2021]() analysing splenic B cells with Il23 knockout which biases B cells away from switches towards the IgG2b isotype.
* [**Comparing cell state transitions inferred using RNA velocity, CSR and SHM**](): This vignette showcases the use of CSR and SHM information to infer B cell state transitions, and compares them to RNA velocity. A subset of tonsillar B cell scRNA-seq data from [King et al. Sci Immunol 2021]() was used as example.
