#===============================================================================
#   Parse phenotypes
#-------------------------------------------------------------------------------

if(!"GEOquery" %in% rownames(installed.packages())){
    source("http://bioconductor.org/biocLite.R")
    biocLite("GEOquery")
}
library(GEOquery)

gse <- getGEO("GSE49031", destdir="data")

all.pheno <- phenoData(gse[[1]])@data
ll <- sapply(all.pheno, function(x) length(levels(x)))
all.meta <- unlist(sapply(all.pheno[1, ll == 1], as.character))
all.pheno <- as.data.frame(lapply(all.pheno[ll > 1], as.character), stringsAsFactors=FALSE)

char.fields <- c(disease.state = "characteristics_ch1",
    immunophenotype = "characteristics_ch1.1",
    subtype = "characteristics_ch1.2",
    cell.type = "characteristics_ch1.3")

all.pheno <- data.frame(
    id = sub(" .*$", "", all.pheno$description),
    geo.accession = all.pheno$geo_accession,
    title = all.pheno$title,
    lapply(char.fields, function(x) factor(sub("^.*: ", "", all.pheno[[x]]))),
    tissue = factor(all.pheno$source_name_ch1),
    stringsAsFactors=FALSE)
idx <- sapply(all.pheno, function(x) "NA" %in% levels(x))
all.pheno[idx] <- lapply(all.pheno[idx], function(x){
    factor(as.character(x), levels=setdiff(levels(x), "NA"))
})
comment(all.pheno) <- all.meta

save(all.pheno, file="data/phenotypes.Rdata")

