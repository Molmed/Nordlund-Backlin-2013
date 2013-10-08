#===============================================================================
#   Download and process gene expression data
#-------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)
require(GEOquery)
require(plyr)

if(!exists("all.pheno")) load("data/phenotypes.Rdata")


#===============================================================================
#   Affymetrix gene expression
#-------------------------------------------------------------------------------

in.file <- "data/GSE47051_series_matrix.txt.gz"
if(!file.exists(in.file))
    download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE47nnn/GSE47051/matrix/GSE47051_series_matrix.txt.gz", in.file)
affy <- getGEO(filename=in.file, destdir="data")


#-----------------------o
#   Expression values

affy.data <- get("exprs", assayData(affy))


#-----------------------o
#   Probe annotations

affy.annot <- with(featureData(affy)@data, data.frame(
    Probe.Set.ID = as.character(ID),
    Gene.Symbol = as.character(get("Gene Symbol")),
    stringsAsFactors = FALSE))


#-----------------------o
#   Phenotypes

affy.pheno <- phenoData(affy)@data
ll <- sapply(affy.pheno, function(x) length(levels(x)))
affy.meta <- unlist(sapply(affy.pheno[1, ll == 1], as.character))

affy.pheno <- data.frame(
    id = sub("^.*(ALL_) ?(\\d+)$", "\\1\\2", affy.pheno$title),
    geo.accession = as.character(affy.pheno$geo_accession),
    title = as.character(affy.pheno$title),
    subtype = factor(sub("^.*: ", "", affy.pheno$characteristics_ch1.2)),
    tissue = affy.pheno$source_name_ch1,
    stringsAsFactors=FALSE)
idx <- sapply(affy.pheno, function(x) "NA" %in% levels(x))
affy.pheno[idx] <- lapply(affy.pheno[idx], function(x){
    factor(as.character(x), levels=setdiff(levels(x), "NA"))
})
comment(affy.pheno) <- affy.meta


#===============================================================================
#   Digital gene expression
#-------------------------------------------------------------------------------

in.file <- "data/GSE26530_DGE_Sense_GeneExpression_TPM.txt.gz"
if(!file.exists(in.file))
    download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE26nnn/GSE26530/suppl/GSE26530_DGE_Sense_GeneExpression_TPM.txt.gz", in.file)
dge <- read.csv(in.file, sep="\t")
dge.geo <- getGEO("GSE26530", destdir="data")


#-----------------------o
#   Expression values

dge.data <- as.matrix(dge[-1])


#-----------------------o
#   Probe annotations

dge.annot <- data.frame(do.call(rbind, strsplit(sub(")$", "", dge$gene.id), "\\(")),
                        stringsAsFactors=FALSE)
names(dge.annot) <- c("transcript.id", "gene.name")


#-----------------------o
#   Phenotypes

# The `dge2met.csv` file will soon be added to the GEO sample series
dge.key <- read.csv("data/dge2met.csv", sep="\t", header=TRUE)
dg <- lapply(dge.geo, phenoData)
dge.pheno <- with(rbind.fill(dg[[1]]@data, dg[[2]]@data), data.frame(
    title = as.character(title),
    geo.accession = as.character(geo_accession),
    stringsAsFactors=FALSE))
dge.pheno <- rev(merge(dge.pheno, dge.key, by.x="title", by.y="dge.id", all.x=TRUE))
names(dge.pheno)[1] <- "id"


#===============================================================================
#   Save
#-------------------------------------------------------------------------------

save(affy.data, affy.annot, affy.pheno, dge.data, dge.annot, dge.pheno, file="data/geneexpression.Rdata")

