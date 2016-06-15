mac_genes <- read.csv("../Cluster/data/Genesets/Macosco_Genesets.csv",sep=";",header = TRUE)

for(x in colnames(mac_genes))
  assign(paste("genes_",x,sep = ""),lapply(mac_genes[[x]][!mac_genes[[x]]==""],function(x) gsub("[[:space:]]","",x)))

