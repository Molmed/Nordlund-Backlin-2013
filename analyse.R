#===============================================================================
#   Setup
#
#   By default this scripts use parallelization through the `foreach`, `SNOW`
#   and `doSNOW` packages. If you cannot run these you will need to manually
#   convert `analyses.R` to single thread mode.
#
    number.of.cores <- 8
#
#-------------------------------------------------------------------------------

tryCatch({
    setwd("data")
    if(!exists("met.data"))  load("methylation.Rdata")
    if(!exists("met.annot")) load("annotations.Rdata")
    if(!exists("met.pheno")) load("phenotypes.Rdata")
    if(!exists("affy.data")) load("geneexpression.Rdata")
    setwd("../results")

    require(doSNOW)
    cl <- makeCluster(number.of.cores)
    registerDoSNOW(cl)
}, error=function(...)
    cat("Could not setup analysis environment. Did you run `setup.R` first?\n"))

site.idx <- with(met.annot, CHR %in% 1:22 & !snp.hit & !bwa.multi.hit)
met.data <- met.data[site.idx,]
met.annot <- met.annot[site.idx,]

types <- table(met.pheno$subtype)
types <- setdiff(names(types)[types >= 10], c("non-recurrent", "undefined"))


#===============================================================================
#   Differentially methylated sites
#-------------------------------------------------------------------------------

ref.idx <- lapply(c(BCP="B", T="T"), function(x)
    grepl(paste0("^(\\d{3}_", x, "|cd34_e1|Constitutional_\\d+)$"), met.pheno$id))
ref.low.var <- lapply(ref.idx, function(x)
    parRapply(cl, met.data[,x], sd, na.rm=TRUE) < .1)

dmc <- pval <- dbeta <- structure(vector("list", length(types)), names=types)
for(my.type in types){
    cat(my.type, "\n")
    my.idx <- met.pheno$subtype %in% my.type
    my.ref <- if(my.type == "T-ALL") "T" else "BCP"
    y <- factor(rep(1:2, c(sum(ref.idx[[my.ref]]), sum(my.idx))),
                labels=c("reference", my.type))
    clusterExport(cl, "y")
    pval[[my.type]] <- parRapply(cl,
        cbind(met.data[,ref.idx[[my.ref]]], met.data[,my.idx]),
        function(x) do.call(wilcox.test, unname(split(x, y)))$p.value)
    dbeta[[my.type]] <- parRapply(cl, met.data[,my.idx], mean, na.rm=TRUE) -
                        parRapply(cl, met.data[,ref.idx[[my.ref]]], mean, na.rm=TRUE)
    dmc[[my.type]] <- p.adjust(pval[[my.type]], "fdr") < .01 &
                      abs(dbeta[[my.type]]) > .2 &
                      ref.low.var[[my.ref]]
}


# Constitutive DMCs
dmc$constitutive <- Reduce("&", dmc[1:8])

# Subtype specific DMCs
dmc[types] <- lapply(dmc[types], "&", Reduce("+", dmc[types]) < 2)

# Relapse DMCs
relapse.ind <- which(met.pheno$sample.type %in% "relapse 1")
diagnosis.ind <- match(sub("r1$", "", met.pheno$id[relapse.ind]),
                       met.pheno$id)
y <- gl(2, length(diagnosis.ind), labels=c("diagnosis", "relapse"))
pval$relapse <- apply(met.data[,c(diagnosis.ind, relapse.ind)], 1,
      function(x) do.call(wilcox.test, c(unname(split(x, y)), list(paired=TRUE)))$p.value)
dbeta$relapse <- apply(met.data[,relapse.ind], 1, mean, na.rm=TRUE) -
                 apply(met.data[,diagnosis.ind], 1, mean, na.rm=TRUE)
dmc$relapse <- pval$relapse < .01 &
               abs(dbeta$relapse) > .2

# Wrap up               
dmg <- lapply(dmc, function(x)
    sort(unique(unlist(strsplit(met.annot$UCSC_REFGENE_NAME[x], ";")))))
sites <- met.annot$TargetID[site.idx]
save(site.idx, sites, dmc, pval, dbeta, dmg, file="dmc.Rdata")


#===============================================================================
#   Correlation with gene expression
#-------------------------------------------------------------------------------

affy.data[affy.data < 3] <- 3
dge.data <- log2(dge.data)
dge.data[dge.data == -Inf] <- 0

sample.ind <- data.frame(id=met.pheno$id,
    met = 1:nrow(met.pheno),
    affy = match(met.pheno$id, affy.pheno$id),
    dge = match(met.pheno$id, dge.pheno$id),
    subtype=as.character(met.pheno$subtype),
    stringsAsFactors=FALSE)
sample.ind$subtype[grepl("^normal", met.pheno$disease.state)] <- "reference"
sample.ind <- subset(sample.ind, (!is.na(affy) | !is.na(dge)) &
                     !subtype %in% c("undefined", NA))

# Reduce dataset, since it will be exported to all workers
site.idx <- Reduce("|", dmc)
dmc <- lapply(dmc, "[", site.idx)
met.data <- met.data[site.idx,]
met.annot <- met.annot[site.idx,]

# Gene to CpG-site mappings
attach(met.annot)
met.genes <- sort(unique(unlist(dmg)))
clusterExport(cl, "UCSC_REFGENE_NAME")
gene2met <- foreach(g = met.genes, .options.mpi=list(chunkSize=20)) %dopar% {
    grep(sprintf("(^|;)%s(;|$)", g), UCSC_REFGENE_NAME)
}
names(gene2met) <- met.genes
met2x <- list(
    dge = foreach(g=met.genes, m=gene2met, .combine=rbind) %dopar% {
        d <- which(dge.annot$gene.name %in% g)
        if(length(d) == 0) return(NULL)
        data.frame(gene=g, expand.grid(met=m, dge=d))
    },
    affy = foreach(g=met.genes, m=gene2met, .combine=rbind) %dopar% {
        a <- grep(paste0("(^| )", g, "( |$)"), affy.annot$Gene.Symbol)
        if(length(a) == 0) return(NULL)
        data.frame(gene=g, expand.grid(met=m, affy=a))
    })
#met2x$dge <- met2x$dge[dge.annot$sense[met2x$dge$dge] %in% "sense",]

types <- table(sample.ind$subtype)
types <- c("all", setdiff(names(types[types > 2]), c("non-recurrent", "reference")))
cor.table <- vector("list", length(types))
names(cor.table) <- types

save.workspace <- function()
    save(sample.ind, site.idx, dmc, dmg, met2x, cor.table,
         envir=globalenv(), file="geneexp.Rdata")
save.workspace()



for(my.type in types){
    cat("Working on", my.type, "\n")

    expr.type <- if(my.type == "all") "dge" else "affy"

    feat <- met2x[[expr.type]][met2x[[expr.type]]$met %in%
        which(dmc[[if(my.type == "all") "constitutive" else my.type]]),]
    feat <- data.frame(feat,
        TargetID = TargetID[feat$met],
        CHR = CHR[feat$met],
        MAPINFO = MAPINFO[feat$met],
        island = RELATION_TO_UCSC_CPG_ISLAND[feat$met])
    feat$region <- mapply(function(g, gg, rr)
        paste(unique(strsplit(rr, ";")[[1]][strsplit(gg, ";")[[1]] %in% g]), collapse=";"),
        feat$gene, UCSC_REFGENE_NAME[feat$met], UCSC_REFGENE_GROUP[feat$met])
    if(expr.type %in% "affy")
        feat$Probe.Set.ID <- affy.annot$Probe.Set.ID[feat$affy]

    if(my.type == "all"){
        group.ind <- !is.na(sample.ind$dge) & !sample.ind$subtype %in% "reference"
        other.ind <- !is.na(sample.ind$dge) & sample.ind$subtype %in% "reference"
    } else {
        group.ind <- !is.na(sample.ind$affy) & sample.ind$subtype %in% my.type
        other.ind <- !is.na(sample.ind$affy) & !sample.ind$subtype %in% my.type
    }
    feat$met.group <- apply(met.data[feat$met, sample.ind$met[group.ind]], 1, mean, na.rm=TRUE)
    feat$met.other <- apply(met.data[feat$met, sample.ind$met[other.ind]], 1, mean, na.rm=TRUE)

    my.expr <- if(expr.type == "dge") dge.data else affy.data
    # Note that this does not copy the data as long as long as we don't modify `my.expr`
    feat$expr.group <- apply(my.expr[feat[[expr.type]], sample.ind[[expr.type]][group.ind]], 1, mean, na.rm=TRUE)
    feat$expr.other <- apply(my.expr[feat[[expr.type]], sample.ind[[expr.type]][other.ind]], 1, mean, na.rm=TRUE)
    feat$log2fold.change <- feat$expr.group - feat$expr.other


    #   ,---.                         |           |
    #   |     ,---. ,---. ,---. ,---. |     ,---. |---  ,---.
    #   |     |   | |     |     |---' |     ,---| |     |---'
    #   `---' `---' `     `     `---' `---' `---^ `---' `---'

    my.met <- data.frame(t(met.data[feat$met, c(sample.ind$met[group.ind], sample.ind$met[other.ind])]))
    my.expr <- data.frame(t(my.expr[feat[[expr.type]], c(sample.ind[[expr.type]][group.ind], sample.ind[[expr.type]][other.ind])]))
    feat$cor.r <- mapply(function(x,y) cor(x, y, use="pair"), my.met, my.expr)


    # Permutation test
    clusterExport(cl, c("my.met", "my.expr"))
    perm.cor <- foreach(icount(10000), .options.mpi=list(chunkSize=20),
                        .inorder=FALSE, .combine=cbind) %dopar% {
        mapply(function(x,y) cor(x, y, use="pair"), my.met[sample(nrow(my.met)),], my.expr)
    }

    feat$cor.p <- apply(abs(cbind(feat$cor.r, perm.cor)), 1, function(x){
        if(is.na(x[1])) return(NA)
        mean(x[1] <= x)
    })
    feat$cor.p.fdr <- p.adjust(feat$cor.p, "fdr")
    cor.table[[my.type]] <- feat
    save.workspace()
}




#===============================================================================
#   Principal component analysis
#-------------------------------------------------------------------------------

##' Regular kNN imputation
##'
##' @param x Dataset.
##' @param k Number of nearest neighbors to use.
##' @param distmat Distance matrix.
##' @return An imputed matrix.
##' @author Christofer \enc{Bäcklin}{Backlin}
##' @export
impute.knn <- function(x, k=.05, distmat){
    if(!is.matrix(x))
        stop("kNN does not work on data with mixed featured types. Therefore as a precausion kNN imputation only accept data in matrix form.")

    if(k < 1) k <- max(1, round(.05*nrow(x)))
    if(k > nrow(x)-1) stop("k is larger than the maximal number of neighbors.")
    if(!is.matrix(distmat)) distmat <- as.matrix(distmat)
    if(any(nrow(x) != dim(distmat)))
        stop("Distance matrix does not match dataset.")

    na.ind <- which(is.na(unname(x)), arr.ind=TRUE)
        # Duplicate names may cause problems otherwise

    NN <- apply(distmat, 1, function(z) order(z))
    fills <- apply(na.ind, 1, function(i){
        mean(na.exclude(x[NN[-1, i[1]], i[2]])[1:k])
    })
    x[na.ind] <- fills
    x
}

met.data <- t(met.data)
all.dist <- dist(met.data)
save(all.dist, file="all_dist.Rdata")

met.data <- impute.knn(met.data, distmat=all.dist)
all.pca <- prcomp(met.data)

plot(all.pca$x[,1], all.pca$x[,2])


