
library(devtools)
library(tensorflow)
library(reticulate)
use_virtualenv(Sys.getenv('ENV'))
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1) 

##source this before sourcing("RunSims.R")
