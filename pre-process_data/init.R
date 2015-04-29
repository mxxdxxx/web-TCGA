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
