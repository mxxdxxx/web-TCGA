require(shiny)
require(shinyIncubator)
require(ggplot2)
require(data.table)
require(markdown)

project.Home="../"

#setwd("~/Documents/BMZ/shinytcga/webapp/")

source("Functions/Helper.R")
sourceDir("Functions/")
sourceDir("tcga.barcode.lib/")

#here we load data
load("r-objects/variant.Data.Table.RData") #snp data
load("r-objects/methylation.Data.Table.RData") # methylation data
#load("r-objects/methylation.Data.Table.RDS.RData") # methylation data
load("r-objects/pos_2_gene.RData") # methylation pos to Gene Symbol
load("r-objects/expression.Data.Table.RData") # Expression Data
load("r-objects/copynumber.Data.Table.RData") # CNV Data

if(system("uname -s",intern=T)=="Darwin"){
  TSS=na.omit(read.table("../Data/tissueSourceSite.csv", sep=",",
                         header=T, stringsAsFactors=F))
} else{
  TSS=na.omit(read.table("/media/daidalos/data/shiny_tcga/tissueSourceSite.csv",
                         sep=",", header=T, stringsAsFactors=F))
}

