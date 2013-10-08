#===============================================================================
#   Download and process infinium methylation array probe annotations
#-------------------------------------------------------------------------------

in.file <- "data/SupplementalTable1.txt.gz"
if(!file.exists(in.file))
    download.file("http://genomebiology.com/imedia/6808530410895336/supp1.gz", "data/SupplementalTable1.txt.gz")

met.annot <- read.csv(in.file,
    header=TRUE, sep="\t", nrows=485577,
    colClasses=c("character", "factor", "integer", "character",
                 "integer", "character", "character", "factor",
                 "character", "character", "character", "character",
                 "factor", "character", "character",
                 # probe filtering
                 "integer", "integer", "integer",
                 # Histone marks
                 "integer", "integer", "integer", "integer", 
                 "integer", "integer", "integer", "integer",
                 # DMCS
                 "integer", "integer", "integer", "integer", 
                 "integer", "integer", "integer", "integer", 
                 "integer", "integer", "integer", "integer"))
# Make the chromosome name a factor manually to get the order right
met.annot$CHR <- factor(met.annot$CHR, levels=c(1:22, c("X", "Y")))
met.annot$CHROMOSOME_36 <- factor(met.annot$CHROMOSOME_36,
    levels=c(levels(met.annot$CHR), "MULTI"))
met.annot[16:38] <- lapply(met.annot[16:38], as.logical)

save(met.annot, file="data/annotations.Rdata")

