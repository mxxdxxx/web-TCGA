project_home="./"
load(paste(project_home,"webapp/r-objects/pat.exp.DT.RData",sep=""))
my.genes=c("MED12","PTEN", "RUNX1", "TP53", "MED12L", "MED13", "MED13L")

types=getCancerTypes(colnames(exp.tumor.mat))

exp.mat=matrix(NA, nrow=length(my.genes), ncol=length(types), dimnames=list(my.genes,types) )
for(i in rownames(exp.mat)){
  for(j in colnames(exp.mat)){
    
    exp.mat[i,j]=res
  }
}