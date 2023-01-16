#' Prepare conda environment
#'
#' @description
#' `prepare_sciCSR` prepares a conda environment with all the python packages required to run sciCSR. You only need to run this function once, before the first time you use the `sciCSR` package.
#'
#' @details
#' `prepare_sciCSR` will set up a conda environment. If Anaconda/Miniconda has not been set up yet in the system, the function will attempt to set up Miniconda using `reticulate::install_miniconda()`.
#' The function will automatically locate Anaconda/Miniconda in the system and create a conda environment containing all dependencies in python which sciCSR uses for its functionalities.
#'
#' @import reticulate
#' @export prepare_sciCSR
prepare_sciCSR <- function()
{
  detect_conda <- try( conda_binary(), silent = TRUE )
  if( inherits(detect_conda, 'try-error') ){
    # set up miniconda
    install_miniconda()
  }
  conda_create(envname = 'scicsr', python_version = '3.9')
  conda_install(envname = 'scicsr', packages = 'scanpy', channel = 'conda-forge')
  conda_install(envname = 'scicsr', packages = 'h5py', pip = TRUE, pip_options = "--force-reinstall")
  conda_install(envname = 'scicsr', packages = 'scvelo', pip = TRUE, pip_options = "-U")
  conda_install(envname = 'scicsr', packages = 'jupyter', channel = 'conda-forge')
  conda_install(envname = 'scicsr', packages = c('python-igraph', 'louvain'),
                pip = TRUE)
  conda_install(envname = 'scicsr', packages = 'cellrank',
                channel = c('conda-forge', 'bioconda'), pip = TRUE)
}