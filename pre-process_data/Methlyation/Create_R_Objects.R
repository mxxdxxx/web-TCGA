# Copy init to run as R script
require(reshape2)
require(markdown)
require(data.table)
#project.Home = "~/Documents/BMZ/shinytcga/"
project.Home = "~/shinytcga/"

TSS = fread(paste(project.Home, "Data/tissueSourceSite.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
TSS = na.omit(TSS)

setnames(TSS, c("TSS.Code","Source.Site","Study.Name", "BCR"))
setkey(TSS, "TSS.Code")

setwd(project.Home)
source("webapp/tcga.barcode.lib/Barcode_Functions.R")




# For all Patients and Genes
# used for "Assembled Normal Probe" Experiment
dir.Content = list.files( paste(project.Home, "Data/Methylation", sep=""),
                          recursive=T)
methylation.Files = dir.Content[grep("*data.txt",dir.Content)]
methylation.Files = paste(project.Home, "Data/Methylation/",
                          methylation.Files, sep="")

cg.ID = fread(methylation.Files[[1]], header=T, skip=1, select=1,
              stringsAsFactors=F)
cg.ID = unlist(cg.ID)

methylation.List.Raw = list()
for(x in methylation.Files){
  for (i in 1:10) gc(reset=T)
  file.Connection = file(x)
  header = readLines(file.Connection, n=2)
  close(file.Connection)
  
  tcga.Barcodes = strsplit(header, "\n")[[1]]
  tcga.Barcodes = strsplit(tcga.Barcodes, "\t")[[1]]
  tcga.Barcodes = tcga.Barcodes[grepl("TCGA", tcga.Barcodes)]
  tcga.Barcodes = unique(tcga.Barcodes)
  length(tcga.Barcodes)
  
  beta.Values.Position = strsplit(header, "\n")[[2]]
  beta.Values.Position = strsplit(beta.Values.Position, "\t")[[1]]
  beta.Values.Position = grep("beta", beta.Values.Position, ignore.case=T)
  
  methylation.File = fread(x, header=T, stringsAsFactors=F,
                           select=beta.Values.Position, skip=1)
  setnames(methylation.File, tcga.Barcodes)
  
  methylation.File = cbind(rn = cg.ID, methylation.File)
  
  #methylation.File = methylation.File[, list(variable = names(.SD),
  #                                           value = unlist(.SD, use.names = F)),
  #                                    by = rn]
  methylation.File = melt(methylation.File)
  setnames(methylation.File, c("rn","patient","beta"))
  setkey(methylation.File, "rn")
  
  methylation.List.Raw[[x]] = methylation.File
  rm(methylation.File)
  for (i in 1:10) gc()
}

# If there is more than one file per entity, merge them -------------------

methylation.List = list()
for(x in methylation.List.Raw){ 
  entity = unique(unlist(getCancerTypes(unique(x$patient))))
  if(is.null(methylation.List[[entity]])){
    methylation.List[[entity]] = rbind(methylation.List[[entity]], x)
  }else{
    methylation.List[[entity]] = x
  }
}

rm(methylation.List.Raw)
save(methylation.List,
     file=paste(project.Home,"webapp/r-objects/methylation.Data.Table.RData" ,sep=""),
     compress="gzip",
     compression_level=4)