###
# We read all gene expression files
# then all meta files
# afterwards we merge them, to get the
# TCGA Barcode-File relation
###

all.uris = list.files(paste(project.Home, "Data/Expression/", sep=""), recursive=T)
exp.Uris = all.uris[grep("mRNAseq_RSEM_normalized_log2.txt", all.uris)]


expression.List = lapply(exp.Uris, function(file){
  exp.File = as.matrix(fread(paste(project.Home, "Data/Expression/" , file, sep=""),
                               header=T, stringsAsFactors=F))
  rownames(exp.File) = exp.File[, 1]
  rownames(exp.File) = gsub("_calculated", "", rownames(exp.File))
  exp.File = exp.File[, -1]
  class(exp.File) = "numeric"
  cols.Keep = (substr(colnames(exp.File), 14, 1) != 1)
  rows.Keep = !grepl("\\?", rownames(exp.File))
  exp.File = exp.File[rows.Keep, cols.Keep]
  new.Rownames = gsub("\\|[0-9]*", "", rownames(exp.File))
  rownames(exp.File) = new.Rownames
  exp.File
})

all.Genes = Reduce(intersect, sapply(expression.List, rownames))

expression.List = lapply(expression.List, function(exp.File){
  keep = which(unique(rownames(exp.File)) %in% all.Genes)
  exp.File = exp.File[keep, ]
})

expression.Table = do.call(cbind, expression.List)
expression.Table = as.data.table(expression.Table, keep.rownames=T)

setkey(expression.Table, "rn")

save(expression.Table, file=paste(project.Home,
                                  "webapp/r-objects/expression.Data.Table.RData",
                                  sep=""),
     compress = "gzip", compression_level = 4)
