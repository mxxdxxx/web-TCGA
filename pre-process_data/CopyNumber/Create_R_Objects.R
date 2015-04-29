dir.Content = list.files( paste(project.Home, "Data/Copynumber", sep=""),
                          recursive=T)
copynumber.Files = dir.Content[grep("all_thresholded.by_genes.txt",dir.Content)]
copynumber.Files = paste(project.Home, "Data/Copynumber/",
                          copynumber.Files, sep="")

gene.Locus.Cyto = list()
copynumber.List = lapply(copynumber.Files, function(x){
  cv.File = fread(x, header=T)
  
  to.Remove = c("Gene Symbol", "Locus ID", "Cytoband")
  gene.Locus.Cyto[[x]] <<- subset(cv.File, T, to.Remove)
  cv.File = subset(cv.File, T, which(!colnames(cv.File) %in% to.Remove))
  cv.File
  })

gene.Locus.Cyto = do.call(rbind, gene.Locus.Cyto)
gene.Locus.Cyto = gene.Locus.Cyto[!duplicated(gene.Locus.Cyto)]
copynumber.Table = do.call(cbind, copynumber.List)
copynumber.Table = cbind(unlist(subset(fread(copynumber.Files[1]),T,1)),
                         copynumber.Table)

colnames(copynumber.Table)[1] = "rn"
setkey(copynumber.Table, "rn")

save(copynumber.Table, gene.Locus.Cyto,
     file=paste(project.Home, "webapp/r-objects/copynumber.Data.Table.RData",
                sep = ""))
