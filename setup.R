#===============================================================================
#
#   This script is the overall pipeline for the analysis code. It downloads all
#   necessary data, sets up the environment and performs the analyses.
#
#   By default only the processed beta values will be downloaded. To also
#   download the unnormalized signal intensities set this variable to TRUE
#
    download.intensities <- FALSE
#
#   By default the analysis scripts use parallelization through the `foreach`,
#   `SNOW` and `doSNOW` packages. If you cannot run these you will need to
#   manually convert `analyses.R` to single thread mode.
#
    number.of.cores <- 8
#
#-------------------------------------------------------------------------------

dir.create("data", showWarnings=FALSE)
dir.create("results", showWarnings=FALSE)

if(file.exists("data/phenotypes.Rdata")){
    load("data/phenotypes.Rdata")
} else {
    source("process_phenotypes.R")
}

if(file.exists("data/annotations.Rdata")){
    load("data/annotations.Rdata")
} else {
    source("process_annotations.R")
}

if(file.exists("data/geneexpression.Rdata")){
    load("data/geneexpression.Rdata")
} else {
    source("process_geneexpression.R")
}

# TODO!!
if(download.intensities && !file.exists("data/intensities.Rdata")){
    source("process_intensities.R")
}

if(file.exists("data/methylation.Rdata")){
    load("data/methylation.Rdata")
} else {
    source("process_methylation.R")
}

source("analyse.R")

