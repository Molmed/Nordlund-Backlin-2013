#===============================================================================
#   Download and procuss infinium methylation array intensities
#
#   This data consist of back ground corrected signal intensities exported from
#   Genome Studio (Illumina inc.).
# #
# #   Line by line instead of using read.table to be more memory efficient.
# #   It is slow, but managable, since it only needs to be done once.
#-------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)
if(!exists("met.annot"))
    load("data/annotations.Rdata")

in.file <- "data/intensities.txt.gz"
if(!file.exists(in.file))
    download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49031/suppl/GSE49031_methylated_unmethylated_signal_intensities.txt.gz", in.file)

# The following section took us about 50 minutes to run, so be patient
columns <- unlist(read.csv(in.file, nrow=1L, sep="\t", header=FALSE))
met.U <- read.csv(in.file, nrow=485577L, header=TRUE, sep="\t",
    colClasses=c(ID_REF="character", Unmethylated="numeric", Methylated="NULL",
                 Detection="NULL")[sub("^\\w+ (\\w+) \\w+$", "\\1", columns)])
if(any(met.U$ID_REF != met.annot$TargetID))
    stop("Probes do not match between annotations and data")
met.U <- as.matrix(met.U[-1])
met.M <- as.matrix(read.csv(in.file, nrow=485577L, header=TRUE, sep="\t",
    colClasses=c(ID_REF="NULL", Unmethylated="NULL", Methylated="numeric",
                 Detection="NULL")[sub("^\\w+ (\\w+) \\w+$", "\\1", columns)]))
met.pval <- as.matrix(read.csv(in.file, nrow=485577L, header=TRUE, sep="\t",
    colClasses=c(ID_REF="NULL", Unmethylated="NULL", Methylated="NULL",
                 Detection="numeric")[sub("^\\w+ (\\w+) \\w+$", "\\1", columns)]))

save(met.U, met.M, met.pval, file="data/intensities.Rdata", compress=FALSE)


#===============================================================================
#   Convert to beta-values and perform peak based correction
#   (Dedeurwaerder et al 2011).
#
#   This is how the data was prepared for the paper.
#-------------------------------------------------------------------------------

# Convert to beta values
met.data <- met.M / (met.U + met.M + 100)
rm(met.U, met.M)
met.data[met.pval > .01] <- NA
rm(met.pval)


# Perform peak-base correction
type <- as.integer(met.annot$INFINIUM_DESIGN_TYPE)
logit <- function(x) log2(x/(1-x))
delogit <- function(x) 2^x/(1+2^x)
get.densities <- function(bw="nrd0"){
    # The type II data in met.data was too large to handle in one go on our computer
    list(density(met.data[type == 1,], bw=bw, na.rm=TRUE),
         density(met.data[which(type == 2)[c(T,F,F)],], bw=bw, na.rm=TRUE),
         density(met.data[which(type == 2)[c(F,T,F)],], bw=bw, na.rm=TRUE),
         density(met.data[which(type == 2)[c(F,F,T)],], bw=bw, na.rm=TRUE))
}
get.peaks <- function(dd, lim){
    if(diff(lim) < 0) stop("Incorrect `lim`")
    sapply(dd, function(d){
        c(with(d, x[x < lim[1]][which.max(y[x < lim[1]])]),
          with(d, x[x > lim[2]][which.max(y[x > lim[2]])]))
    })
}

met.data <- logit(met.data)
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
    met.data[i,] <- met.data[i,] * typeII.scale[1 + (met.data[i,] >= 0)]
}
met.data <- delogit(met.data)

save(met.data, file="data/methylation.Rdata", compress=FALSE)


