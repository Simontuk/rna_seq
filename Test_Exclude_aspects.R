
#Excluding Aspect "GSEA_c2_Rosty..."
test <- (tamr2$cnam$`#PC1# GSEA_c2_ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER`)

gsub("#PC1# ","",names(test))
cc.pattern <- pagoda.show.pathways((gsub("#PC1# ","",names(test))), varinfo, go.env, show.cell.dendrogram = TRUE, cell.clustering = hc, showRowLabels = TRUE)

varinfo.cc <- pagoda.subtract.aspect(varinfo, cc.pattern)
head(varinfo.cc$mat[,1:5])

pwpca <- pagoda.pathway.wPCA(varinfo.cc, go.env, n.components = 1, n.cores = 12)

df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)

tam <- pagoda.top.aspects(pwpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))

tamr <- pagoda.reduce.loading.redundancy(tam, pwpca)

tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)

col.cols <- rbind(groups = cutree(hc, 4))

pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20))
