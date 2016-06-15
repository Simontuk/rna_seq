## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library("scde")
library("psych")
library("biomaRt")
library("knitr")
library("GSEABase")
library("RColorBrewer")
library("gplots")


## ----Data Summary--------------------------------------------------------
load("~/RNA_Seq/For_Simon/cnag_count_all_experiments.RData")
ls()
head(summary(count_list))
kable(head(summary(count_list$MCF10CA_3D_5d_random[,1:3])))


## ----list_Samples, echo = FALSE------------------------------------------
rownames(head(summary(count_list)))
MCFsample <- "MCF10A_3D_5d_random"


## ----Output Table Countlist, echo=FALSE----------------------------------
kable(count_list$MCF10CA_2D_random[1:5,1:4])


## ----Output Table NomrVal, echo=FALSE------------------------------------
kable(deseq2_list$MCF10CA_2D_random[1:5,1:4])


## ----Data preparation----------------------------------------------------
MCF_rawcounts <- count_list[[MCFsample]]
colnames(MCF_rawcounts)


## ----CleanCounts---------------------------------------------------------
ccounts <- clean.counts(MCF_rawcounts, min.lib.size = 1000, min.reads = 1, min.detected = 1)
dim(MCF_rawcounts)
dim(ccounts)
rownames(ccounts) <- gsub("\\.\\d+", "\\1", rownames(ccounts)) # Removing Ensembl_ID decimal (xxx.4)
ccounts <- ccounts[,1:69] # Removing the "Bulk Experiment"


## ----Fitting Error Model, eval = FALSE-----------------------------------
#o.ifm <- scde.error.models(counts = ccounts, groups = NULL, n.cores = 8, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
#save(o.ifm, file=("~/RNA_Seq/For_Simon/errormodelfit.rdata"))
## 

## ----Load Error Mode-----------------------------------------------------
load(file=("~/RNA_Seq/For_Simon/errormodelfit.rdata"))


## ----Filtering cells-----------------------------------------------------
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)

o.ifm <- o.ifm[valid.cells, ]


## ----Estimate gene expression, eval=FALSE,echo=FALSE---------------------
## o.prior <- scde.expression.prior(models = o.ifm, counts = ccounts, length.out = 400, show.plot = FALSE)
## 

## ----DiffExp Testing, eval=FALSE,echo=FALSE------------------------------
## ediff <- scde.expression.difference(o.ifm, ccounts, o.prior, n.randomizations = 100, n.cores = 8, verbose = 1)
## 

## ----KnearestNeighbor Error Models, eval = FALSE-------------------------
## # not calculated due to high complexity, but sample plots are included
ko.ifm <- knn.error.models(counts = ccounts, k = ncol(ccounts)/4, n.cores = 7, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)
## 

## ----Normalize Variance, eval = FALSE------------------------------------
## # not calculated due to high complexity, resulting plot is included
## varinfo <- pagoda.varnorm(ko.ifm, counts = ccounts, trim = 3/ncol(ccounts), max.adj.var = 5, n.cores = 8, plot = TRUE)
## save(list = c("ko.ifm","varinfo"), file=("~/RNA_Seq/For_Simon/KNNerrormodelfit.rdata"))
## 

## ----Load_Var/knn_Data, echo=FALSE---------------------------------------
load(file=("~/RNA_Seq/For_Simon/KNNerrormodelfit.rdata"))

## ----top overdispersed genes---------------------------------------------
# list top overdispersed genes
sort(varinfo$arv, decreasing = TRUE)[1:10]


## ----Controlling for sequencing depth------------------------------------
varinfo <- pagoda.subtract.aspect(varinfo, colSums(ccounts[, rownames(ko.ifm)]>0))


## ----Package org.Hs.eg.db, include=FALSE---------------------------------
library(org.Hs.eg.db)

## ----Creating go-term environment----------------------------------------
# translate gene names to EG_ids
ids <- unlist(lapply(mget(rownames(ccounts), org.Hs.egENSEMBL2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids); names(rids) <- ids

# convert GO lists from ids to gene names
gos.interest <- unique(c(ls(org.Hs.egGO2ALLEGS)[1:100],"GO:0022008","GO:0048699", "GO:0000280", "GO:0007067"))

go.env <- lapply(mget(gos.interest, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x])))


----Other Gene Sets, eval=FALSE-----------------------------------------
# Read in a list of gene sets
genesetsNames <- read.csv("~/RNA_Seq/RNA-Seq-Project/Testgenesets.csv",header = FALSE,sep = ",")
genuri <- asBroadUri(genesetsNames) # Convert to URI
##
genesetsofinterest <- getBroadSets(genuri) # Download gene sets
##
# Create an Environment with the Gene sets and included genes
gois <- sapply(names(genesetsofinterest), function(x) geneIds(genesetsofinterest[[x]]), USE.NAMES = TRUE)
##
go.env <- clean.gos(gois) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment
##
##
#Change rownames to Gene-Symbols
ids <- unlist(lapply(mget(rownames(ccounts), org.Hs.egENSEMBL2EG, ifnotfound = NA), function(x) x[1]))
##
id1 <- na.omit(ids)
##
symbs <- unlist(lapply(mget(id1, org.Hs.egSYMBOL, ifnotfound = NA), function(x) x[1]))
##
rids <- symbs; names(rids) <- names(id1)
nasymbs <- unlist(lapply(rownames(ccounts), function(x) rids[x])) # Switching EnsemblID to GeneSymbol
##
ccountstest <- ccounts
rownames(ccountstest) <- nasymbs
##
ccountstest <- ccountstest[!is.na(rownames(ccountstest)),]
##
##
# Error Model for new Element with gene symbols
##
ko.ifm <- knn.error.models(counts = ccountstest, k = ncol(ccountstest)/4, n.cores = 8, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)
##
# Normalization
##
varinfo <- pagoda.varnorm(ko.ifm, counts = ccountstest, trim = 3/ncol(ccountstest), max.adj.var = 5, n.cores = 8, plot = TRUE)
##
# Controlling for sequencing depth
varinfo <- pagoda.subtract.aspect(varinfo, colSums(ccountstest[, rownames(ko.ifm)]>0))
##
# First principal component analysis
##
pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 8)
##
##
##

----Testing, eval=FALSE-------------------------------------------------
hall_e_m <- getBroadSets(uri = asBroadUri(c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")))
##
GeneSet(hall_e_m,geneIdType=EntrezIdentifier("org.Hs.eg.db"),setName="Sample")
hall_e_m <- mapIdentifiers(hall_e_m, EntrezIdentifier())
##
##
##
gidsHall <- geneIds(hall_e_m)
##
gidsHall <- gidsHall$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
##
Hall_list <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","ensembl_gene_id", "entrezgene"),values=gidsHall, mart=mart)
##
##
showMethods("GeneSetCollection")
##
geneIds(hall_e_m)
##
##
gos.interest.a <- c("GO:0007067","GO:0007049")
##
##
##
##
go.env.a <- lapply(mget(gos.interest.a, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x])))
##
go.env <- go.env.a
##
go.env <- NA
go.env$'Hallmark_Mes_Ench' <- Hall_list$ensembl_gene_id

## ----env clean up--------------------------------------------------------
# clean up and conversion to environment:
go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment


## ----Weighted first PCA for each GO gene set-----------------------------
pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 8)


## ----Statistical significance--------------------------------------------
df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)


## ----Table Overdispersed gene sets---------------------------------------
kable(head(df))


## ----Calculate gene clusters, eval=FALSE---------------------------------
## # This does not work on local machine
## clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 2, n.cores = 4, plot = TRUE)
##
##

## ----Quick Webinterface--------------------------------------------------
tam <- pagoda.top.aspects(pwpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))

hc <- pagoda.cluster.cells(tam, varinfo)

tamr <- pagoda.reduce.loading.redundancy(tam, pwpca)

tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)


col.cols <- rbind(groups = cutree(hc, 3), group2 = cutree(hc,7))
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20))


rownames(tam$xv) <- unlist(lapply(strsplit(rownames(tam$xv),split = "E_"),function(x) x[2]))
hmcol<-brewer.pal(11,"RdBu")
heatmap.2(tam$xv, trace = "none", col = hmcol,colCol =brewer.pal(5,name = "Set1")[cutree(hc,5)])


# compile a browsable app, showing top three clusters with the top color bar
app <- make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, col.cols = col.cols, cell.clustering = hc, title = "NPCs")

save(app, file=("~/RNA_Seq/For_Simon/app.rdata"))


## ----prompt=TRUE---------------------------------------------------------
system("Rscript OpenWebApp.R")
## Use in R in terminal:
# show app in the browser (port 1468)


## ----Scrapyard, eval=FALSE-----------------------------------------------
## ## Playground:
## mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
##
## genes<-c("ENSG00000000003.10","ENSG00000000005.5","ENSG00000000419.8","ENSG00000000457.8","ENSG00000000460.11")
##
##
## genes <- gsub("\\.\\d+", "\\1", genes)
## genes
##
## library(org.Hs.eg.db)
##
## ids <- unlist(lapply(mget(rownames(ccounts), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
## rids <- names(ids); names(rids) <- ids
##
##
## getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene","go_id"),values=genes,mart= mart)
##
## ## Reducing the number of evaluated cells to test `pagoda.gene.clusters` function
## dim(ccounts)
## ccounts1 <- ccounts[,1:20]
##
## ko.ifm1 <- knn.error.models(counts = ccounts1, k = ncol(ccounts)/4, n.cores = 8, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10,save.model.plots = FALSE)
##
## varinfo1 <- pagoda.varnorm(ko.ifm1, counts = ccounts1, trim = 3/ncol(ccounts), max.adj.var = 5, n.cores = 8, plot = TRUE)
##
## varinfo1 <- pagoda.subtract.aspect(varinfo1, colSums(ccounts1[, rownames(ko.ifm1)]>0))
##
## clpca1 <- pagoda.gene.clusters(varinfo1, trim = 7.1/ncol(varinfo1$mat), n.clusters = 10, n.cores = 8, plot = TRUE)
##
##
##
##
##
##

----Reduced Gene Set, eval=FALSE----------------------------------------
#Change rownames to Gene-Symbols
ids <- unlist(lapply(mget(rownames(ccounts), org.Hs.egENSEMBL2EG, ifnotfound = NA), function(x) x[1]))
##
id1 <- na.omit(ids)
##
symbs <- unlist(lapply(mget(id1, org.Hs.egSYMBOL, ifnotfound = NA), function(x) x[1]))
##
rids <- symbs; names(rids) <- names(id1)
nasymbs <- unlist(lapply(rownames(ccounts), function(x) rids[x])) # Switching EnsemblID to GeneSymbol
##
ccountstest <- ccounts
rownames(ccountstest) <- nasymbs
##
ccountstest <- ccountstest[!is.na(rownames(ccountstest)),]
##
##
# Error Model for new Element with gene symbols
##
ko.ifm <- knn.error.models(counts = ccountstest, k = ncol(ccountstest)/4, n.cores = 8, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)
##
# Normalization
##
varinfo <- pagoda.varnorm(ko.ifm, counts = ccountstest, trim = 3/ncol(ccountstest), max.adj.var = 5, n.cores = 8, plot = TRUE)
##
# Controlling for sequencing depth
varinfo <- pagoda.subtract.aspect(varinfo, colSums(ccountstest[, rownames(ko.ifm)]>0))
##
##
# Read in a list of gene sets
#genesetsNames <- read.csv("~/RNA_Seq/RNA-Seq-Project/Testgenesets.csv",header = FALSE,sep = ",")
#genuri <- asBroadUri(genesetsNames) # Convert to URI
##
#genesetsofinterest <- getBroadSets(genuri) # Download gene sets
##
# Create an Environment with the Gene sets and included genes
#gois <- sapply(names(genesetsofinterest), function(x) geneIds(genesetsofinterest[[x]]), USE.NAMES = TRUE)
##
cyclebase.peaktime <- read.csv("../CycleBase_files/human_periodic.tsv",sep = "\t")
##
Periodic_Cyclebase <- cyclebase.peaktime$gene
##
Periodic_Cyclebase <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_gene_id", "entrezgene","hgnc_symbol","ensembl_peptide_id"),values=Periodic_Cyclebase,mart= mart)
cyclebase.peaktime <- cyclebase.peaktime[cyclebase.peaktime$gene %in% Periodic_Cyclebase$ensembl_peptide_id,]
##
Periodic_Cyclebase <- Periodic_Cyclebase[!duplicated(Periodic_Cyclebase$ensembl_peptide_id,),]
cyclebase.peaktime$gensymb <- Periodic_Cyclebase$hgnc_symbol
##
gois <- NULL
gois$cyclebase <- Periodic_Cyclebase$hgnc_symbol
###################
# Read Whitfield Genesets:
genesetsNames <- read.csv("~/RNA_Seq/RNA-Seq-Project/Testgenesets.csv",header = FALSE,sep = ",")
genuri <- asBroadUri(genesetsNames) # Convert to URI
##
genesetsofinterest <- getBroadSets(genuri) # Download gene sets
##
# Create an Environment with the Gene sets and included genes
gois_wf <- sapply(names(genesetsofinterest), function(x) geneIds(genesetsofinterest[[x]]), USE.NAMES = TRUE)
##
gois_wf[["G1&S"]] <- c(gois_wf$WHITFIELD_CELL_CYCLE_G1_S,gois_wf$WHITFIELD_CELL_CYCLE_S)
gois_wf[["G2&M"]] <- c(gois_wf$WHITFIELD_CELL_CYCLE_G2,gois_wf$WHITFIELD_CELL_CYCLE_G2_M)
gois_wf[["M_G1"]] <- gois_wf$WHITFIELD_CELL_CYCLE_M_G1
###################
##
gois[["G1&S"]] <- gois_wf[["G1&S"]]
gois[["G2&M"]] <- gois_wf[["G2&M"]]
gois[["M_G1"]] <- gois_wf[["M_G1"]]
##
go.env <- clean.gos(gois) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment
##
##
##
# First principal component analysis

pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 8)


tam <- pagoda.top.aspects(pwpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))

hc <- pagoda.cluster.cells(tam, varinfo)

tamr <- pagoda.reduce.loading.redundancy(tam, pwpca,plot = TRUE)

tamr2 <- pagoda.reduce.redundancy(tam, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = rownames(tam$xv), labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0.2)


col.cols <- rbind(groups = cutree(hc, 3), group2 = cutree(hc,7))
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20))


rownames(tam$xv)

rownames(tam$xv) <- unlist(lapply(strsplit(rownames(tam$xv),split = " "),function(x) x[2]))
hmcol<-brewer.pal(11,"RdBu")
heatmap.2(tam$xv, trace = "none", col = hmcol,colCol =brewer.pal(5,name = "Set1")[cutree(hc,5)])


## ------------------------------------------------------------------------

cyclebase.peaktime <- cyclebase.peaktime[cyclebase.peaktime$gene %in% Periodic_Cyclebase$ensembl_peptide_id,]

gois$cyclebase <- Periodic_Cyclebase$hgnc_symbol

hmcol<-brewer.pal(11,"RdBu")

#Gene / Cell clustering for cyclebase geneset
i <- 4 # N of bins for Peaktime
countscycle <- (varinfo$mat[rownames(varinfo$mat) %in% gois$cyclebase,])
sidecol <- brewer.pal(i,"Greens")
countsycle.peaktime <- cyclebase.peaktime[cyclebase.peaktime$gensymb %in% rownames(countscycle),]
cutted <- cut(countsycle.peaktime$peaktime,i)
heatmap.2(countscycle, trace = "none", col = hmcol, RowSideColors = sidecol[cutted[-1]],revC = TRUE)

countscycle <- (varinfo$mat[rownames(varinfo$mat) %in% gois$`G1&S`,])
heatmap.2(countscycle, trace = "none", col = hmcol,key.title = "G1&S",colCol =brewer.pal(5,name = "Set1")[cutree(hc,3)])

countscycle <- (varinfo$mat[rownames(varinfo$mat) %in% gois$`G2&M`,])
heatmap.2(countscycle, trace = "none", col = hmcol,key.title = "G2&M", colCol =brewer.pal(5,name = "Set1")[cutree(hc,3)])

countscycle <- (varinfo$mat[rownames(varinfo$mat) %in% gois$M_G1,])
heatmap.2(countscycle, trace = "none", col = hmcol,key.title = "M_G1", colCol =brewer.pal(5,name = "Set1")[cutree(hc,3)])

countscycle <- (varinfo$mat[rownames(varinfo$mat) %in% gois_wf$WHITFIELD_CELL_CYCLE_LITERATURE,])
heatmap.2(countscycle, trace = "none", col = hmcol,colCol =brewer.pal(5,name = "Set1")[cutree(hc,3)])

colnames(varinfo$mat)


