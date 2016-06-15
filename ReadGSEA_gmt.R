fl <- list.files(path = "../GSEA_genesets",full.names = TRUE)

GeneCollections <- NA 
for(x in 1:6)
  GeneCollections[x] <- paste(sep = "","GSEA_",(strsplit(strsplit(fl,split = "/")[[x]][3],split = "\\.")[[1]][1]))

for(x in 1:6)
  assign(GeneCollections[x],getGmt(fl[x]))


gois <- sapply(GeneCollections, function(x) 
                    geneIds(get(x)), USE.NAMES = TRUE)


for(x in names(gois))
  for(i in names(gois[[x]]))
    gois[[paste(x,i,sep="_")]] <- gois[[x]][[i]]

go.env <- clean.gos(gois) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment
