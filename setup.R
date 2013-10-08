#===============================================================================
#
#   This script sets up the environment required for the analysis code.
#   After this script has finised you can run `analyse.R`.
#
#   By default only the processed beta values will be downloaded. To also
#   download the unnormalized signal intensities set this variable to TRUE
#
    download.intensities <- FALSE
#
#-------------------------------------------------------------------------------


#-------------------------------o
#   Make directory structure

dir.create("data", showWarnings=FALSE)
dir.create("results", showWarnings=FALSE)


#-------------------------------o
#   Install required packages

required.pkg <- c("plyr", "doSNOW")
required.bioc.pkg <- c("GEOquery")
installed.pkg <- rownames(installed.packages())

required.pkg <- required.pkg[!required.pkg %in% installed.pkg]
required.bioc.pkg <- required.bioc.pkg[!required.bioc.pkg %in% installed.pkg]

if(!is.null(required.pkg)) install.packages(required.pkg)
if(!is.null(required.bioc.pkg)){
    source("http://bioconductor.org/biocLite.R")
    biocLite(required.bioc.pkg)
}


#-------------------------------o
#   Download and prepare data

if(!file.exists("data/phenotypes.Rdata")){
    source("process_phenotypes.R")
}

if(!file.exists("data/annotations.Rdata")){
    source("process_annotations.R")
}

if(!file.exists("data/geneexpression.Rdata")){
    source("process_geneexpression.R")
}

# TODO!!
if(download.intensities && !file.exists("data/intensities.Rdata")){
    source("process_intensities.R")
}

if(!file.exists("data/methylation.Rdata")){
    source("process_methylation.R")
}

