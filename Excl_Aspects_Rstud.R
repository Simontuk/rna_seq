test <- read.csv("../Cluster/data/AspectsToExclude.csv",sep=";",header=FALSE)
rownames(test) <- test[,1]
test <- test[,-1]


PC1names <- names(tamr2$cnam)


exclAsp <- na.exclude(unlist(test["CA_3D_5d_pooled_invasive",],use.names = FALSE))


tamr2$cnam$`#PC1# GSEA_c2_ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER`[1:5]


nli <- unlist(PC1names[exclAsp])

pws_list <- unlist(sapply(nli, function(x) tamr2$cnam[[x]],USE.NAMES = FALSE))


Excl_pws <- gsub("#PC1# ","",names(pws_list))
Excl_pws <- gsub("1$","",Excl_pws)

cc.pattern <- pagoda.show.pathways(Excl_pws, varinfo, go.env, show.cell.dendrogram = TRUE, cell.clustering = hc, showRowLabels = TRUE)


varinfo.cc <- pagoda.subtract.aspect(varinfo, cc.pattern)

