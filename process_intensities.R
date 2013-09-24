#===============================================================================
#   Parse annotations
#-------------------------------------------------------------------------------

if(!"GEOquery" %in% rownames(installed.packages())){
    source("http://www.bioconductor.org/biocLite.R")
    biocLite("GEOquery")
}
require("GEOquery")


x <- getGEO("GSE49031")

download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49031/suppl/GSE49031_methylated_unmethylated_signal_intensities.txt.gz", "intensities.txt.gz")





#===============================================================================
#   Download and parse infinium methylation array intensities into R
#
#   This data consist of back ground corrected signal intensities exported from
#   Genome Studio (Illumina inc.).
#
#   Line by line instead of using read.table to be more memory efficient.
#   It is slow, but managable, since it only needs to be done once.
#-------------------------------------------------------------------------------

in.file <- "intensities.txt"

con <- file(in.file, "r")
read.con <- function() strsplit(readLines(con, n=1L), "\t")[[1]]
columns <- read.con()
ind <- grep("Unmethylated Signal$", columns)
all.met <- matrix(NA, 485577, length(ind),
    dimnames=list(site=met.annot$TargetID,
        sample=sub(" Unmethylated Signal$", "", columns[ind])))
site.names <- rep(NA, nrow(all.met))

for(i in 1:nrow(all.met)){
    l <- read.con()
    site.names[i] <- l[1]
    a <- as.numeric(l[ind])
    b <- as.numeric(l[ind+1])
    # Remove probles with detection p-value > .01
    all.met[i,] <- ifelse(as.numeric(l[ind+2]) <= 0.01, b/(a+b+100), NA)
}
close(con)
if(any(site.names != met.annot$TargetID))
    stop("Probes do not match between annotations and data")


#===============================================================================
#   Peak based correction (Dedeurwaerder et al 2011)
#-------------------------------------------------------------------------------

type <- as.integer(met.annot$INFINIUM_DESIGN_TYPE)
logit <- function(x) log2(x/(1-x))
delogit <- function(x) 2^x/(1+2^x)
get.densities <- function(bw="nrd0"){
    # The type II data in met.data is too large to be handled in one go
    list(density(all.met[type == 1,], bw=bw, na.rm=TRUE),
         density(all.met[which(type == 2)[c(T,F,F)],], bw=bw, na.rm=TRUE),
         density(all.met[which(type == 2)[c(F,T,F)],], bw=bw, na.rm=TRUE),
         density(all.met[which(type == 2)[c(F,F,T)],], bw=bw, na.rm=TRUE))
}
get.peaks <- function(dd, lim){
    if(diff(lim) < 0) stop("Incorrect `lim`")
    sapply(dd, function(d){
        c(with(d, x[x < lim[1]][which.max(y[x < lim[1]])]),
          with(d, x[x > lim[2]][which.max(y[x > lim[2]])]))
    })
}

all.met <- logit(all.met)
# Since we're only interested in the location of the modes we can use a
# predefined bandwidth, which is faster to compute but coverges more slowly
# towards the true distribution.
den.rescaling <- get.densities(.5)
peaks.rescaling <- get.peaks(den.rescaling, c(-1.5, 1.5))

# Scale M-values and calculate peaks
typeII.scale <- c(lower = peaks.rescaling[1,1] / mean(peaks.rescaling[1,-1]),
                  upper = peaks.rescaling[2,1] / mean(peaks.rescaling[2,-1]))
# A vectorized solution would be faster, but for loops require less memory
for(i in which(type == 2)){
    all.met[i,] <- all.met[i,] * typeII.scale[1 + (all.met[i,] >= 0)]
}
all.met <- delogit(all.met)



#===============================================================================
#   Export the processed data, both as .Rdata and .txt files
#-------------------------------------------------------------------------------

in.file <- "intensities.txt"
out.file <- "intensities_processed.txt"
save(all.met, file=sub("\\.txt", ".Rdata", out.file))

con <- file(in.file, "r")
read.con <- function() strsplit(readLines(con, n=1L), "\t")[[1]]
my.cat <- function(..., append=TRUE){
    cat(..., sep="\t", file=out.file, append=append)
    cat("\n", file=out.file, append=TRUE)
}
columns <- read.con()
ind <- grep("Detection Pval$", columns)

my.cat("ID_REF",
       rbind(sub("Unmethylated Signal$", "Average Beta", columns[ind]),
          sub("Unmethylated Signal$", "Detection Pval", columns[ind])), append=FALSE)
for(i in 1:nrow(all.met)){
    l <- read.con()
    my.cat(l[1], rbind(sprintf("%.7f", all.met[i,]), l[ind]))
}
close(con)

