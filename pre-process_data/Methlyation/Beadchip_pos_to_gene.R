pos.2.gene=fread(paste(project_home,"/Data/humanmethylation450_15017482_v1-2.csv",sep=""),header=T,stringsAsFactors=F)
keep.columns=c("IlmnID","Chromosome_36","Coordinate_36","Probe_SNPs_10","UCSC_RefGene_Name", "UCSC_RefGene_Group", "UCSC_CpG_Islands_Name")
pos.2.gene=pos.2.gene[,keep.columns,with=F]

genes=strsplit(pos.2.gene$UCSC_RefGene_Name,";")
region=strsplit(pos.2.gene$UCSC_RefGene_Group,";")

genes=sapply(genes,function(x) x[1])
region=sapply(region,function(x) x[1])

pos.2.gene$UCSC_RefGene_Name=genes
pos.2.gene$UCSC_RefGene_Group=region

setkey(pos.2.gene, "UCSC_RefGene_Name")
save(pos.2.gene,file=paste(project_home,"webapp/r-objects/pos_2_gene.RData" ,sep=""), compress="gzip", compression_level=4)
