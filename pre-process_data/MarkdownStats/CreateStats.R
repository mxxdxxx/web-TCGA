require(knitr)
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
load("r-objects/expression.Data.Table.RData") # Expression Data
load("r-objects/copynumber.Data.Table.RData") # CNV Data

knit(paste(project.Home, "pre-process_data/MarkdownStats/CreateStats.Rmd", sep=""),
         output=paste(project.Home, "webapp/r-objects/Stats.md", sep=""),quiet=T)