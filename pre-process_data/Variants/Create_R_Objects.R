
variant_path = paste(project.Home, "Data/Variant/", sep="")
all.files = list.files(variant_path,recursive = T)
maf.files = all.files[grepl("maf", all.files)]

mafs = list()
for(i in 1:length(maf.files)){
  f = maf.files[[i]]
  path_to_read = paste(variant_path, f, sep = "")
  maf = read.csv(path_to_read, sep = "\t", header = T)
  maf = subset(maf, T, c(1:14, 16, 17))
  mafs[[i]] = as.data.table(maf)
}

mafs = rbindlist(mafs)
#mafs = subset(mafs, subset = T, select = c(1:14, 16, 17))

CTs = getTSS(unique(mafs$Tumor_Sample_Barcode))
CTs = cbind(CTs,as.character(unique(mafs$Tumor_Sample_Barcode)))
colnames(CTs) = c("Entity","Tumor_Sample_Barcode")
mafs = merge(mafs, CTs, by = "Tumor_Sample_Barcode")
variant.Table = mafs
variant.Table$Entity = as.character(variant.Table$Entity)

setnames(variant.Table,which(grepl("norm", colnames(variant.Table),
                                   ignore.case = T)), "Normal_Sample_Barcode")

variant.Table = cbind(variant.Table, paste(variant.Table$Hugo_Symbol,": ",
                                           variant.Table$Variant_Classification,
                                           sep = "" ))
setnames(variant.Table, ncol(variant.Table), "Combined_Class")

variant.Table$Tumor_Sample_Barcode = as.character(variant.Table$Tumor_Sample_Barcode)
variant.Table$Hugo_Symbol = as.character(variant.Table$Hugo_Symbol)
variant.Table$Chromosome = as.character(variant.Table$Chromosome)
variant.Table$Strand = as.character(variant.Table$Strand)
variant.Table$Entity = as.character(variant.Table$Entity)
variant.Table$Combined_Class = as.character(variant.Table$Combined_Class)

setkey(variant.Table, "Hugo_Symbol")

save(variant.Table,
     file = paste(project.Home, "webapp/r-objects/variant.Data.Table.RData", 
                  sep = ""), compress = "gzip", compression_level = 4)