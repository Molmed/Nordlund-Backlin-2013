#===============================================================================
#   Download and import methylation into R
#
#   This data consist of beta values normalized with peak based correction
#   (Dedeurwaerder et al 2011).
#
#   Beware that you need approximately 20 GB of RAM to run this script
#   efficiently, since the gzipped text file is read into a 6.9 GB data frame
#   (x) and processed into a 3.4 GB numeric matrix as the final result, leaving
#   some computational headroom.
#-------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)
if(!exists("met.annot"))
    load("data/annotations.Rdata")

in.file <- "data/processed.txt.gz"
if(!file.exists(in.file))
    download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49031/suppl/GSE49031_processed.txt.gz", in.file)

columns <- unlist(read.csv(in.file, nrow=1L, sep="\t", header=FALSE))
x <- read.csv(in.file, nrow=485577L, header=TRUE, sep="\t",
    colClasses=c(ID_REF="character", Beta="numeric", Pval="numeric")[sub(".* ", "", columns)])

all.met <- do.call(cbind, x[grep("Beta", names(x))])
all.met.pval <- do.call(cbind, x[grep("Pval", names(x))])
rm(x)
all.met[all.met.pval > .01] <- NA
rm(all.met.pval)
all.met <- all.met[met.annot$analyzed,]

save(all.met, file="data/methylation.Rdata", compress=FALSE)

