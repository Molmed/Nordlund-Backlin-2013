#===============================================================================
#   Setup
#-------------------------------------------------------------------------------

if(!file.exists("results")) dir.create("results")
if(!exists("number.of.cores")) number.of.cores <- 8

library(doSNOW)
cl <- makeCluster(number.of.cores)
registerDoSNOW(cl)

setwd("data")
if(!exists("met.pheno")) load("phenotypes.Rdata")
if(!exists("met.annot")) load("annotations.Rdata")
if(!exists("met.data"))  load("methylation.Rdata")
if(!exists("affy.data")) load("geneexpression.Rdata")
setwd("../results")




site.idx <- with(met.annot, CHR %in% 1:22 & !snp.hit & !bwa.multi.hit)
met.data <- met.data[site.idx,]
met.annot <- met.annot[site.idx,]

patterns <- c(control="(MethDNA(Pos|Neg)_\\d+|WGA_CTRL_\\d+)",
    reference="(\\d{3}_[BT]|cd34_e1)", diagnosis="ALL_\\d+",
    remission="Constitutional_\\d+", `relapse 1`="ALL_\\d+r1", 
    `relapse 2`="ALL_\\d+r2")
met.pheno$sample.type <- factor(NA, levels=names(patterns))
for(i in seq(patterns))
    met.pheno$sample.type[grep(patterns[i], met.pheno$sample.name)] <- names(patterns)[i]

types <- table(met.pheno$subtype)
types <- setdiff(names(types)[types >= 10], "undefined")


#===============================================================================
#   Differentially methylated sites
#-------------------------------------------------------------------------------

ref.idx <- lapply(c(BCP="B", T="T"), function(x)
    grepl(paste0("^(\\d{3}_", x, "|cd34_e1|Constitutional_\\d+)$"), met.pheno$sample.name))
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
diagnosis.ind <- match(sub("r1$", "", met.pheno$sample.name[relapse.ind]),
                       met.pheno$sample.name)
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

library(analyse450k)
load.450k.anot(c("affy", "dge"))

affy.data[affy.data < 3] <- 3
dge.data <- log2(dge.data)
dge.data[dge.data == -Inf] <- 0

# Sample mappings
dge.id <- dge.pheno$id
i <- dge.pheno$sample.type %in% c("healthy B", "healthy T")
dge.id[i] <- sub("([BT])_(\\d+)", "\\2_\\1", dge.id[i])
dge.id[!i] <- as.character(dge.pheno$fmca.id[!i])

sample.ind <- do.call(rbind, lapply(met.pheno$id, function(my.id){
    data.frame(id=my.id,
               met=which(met.pheno$id %in% my.id),
               dge=tail(c(NA, which(dge.id %in% my.id)), 1),
               affy=c(which(affy.pheno$id %in% my.id), NA)[1],
               subtype=met.pheno$subtype[met.pheno$id %in% my.id])
}))
sample.ind$subtype <- as.character(sample.ind$subtype)
sample.ind$subtype[grepl("^normal", met.pheno$accession)] <- "reference"
sample.ind <- sample.ind[(!is.na(sample.ind$dge) | !is.na(sample.ind$affy))
                         & !sample.ind$subtype %in% c("undefined", NA),]

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
met2x$dge <- met2x$dge[dge.annot$sense[met2x$dge$dge] %in% "sense",]

types <- table(sample.ind$subtype)
types <- c("all", names(types[types > 2]))
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


