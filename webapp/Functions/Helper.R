head.df=function(df){
  if(nrow(df) < 5 & ncol(df) < 5) return(df)
  if(nrow(df) < 5 ) return(df[,1:5])
  if(ncol(df) < 5 ) return(df[1:5,])
}

sourceDir = function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

#Split gene names by comma
#remove white spaces
process_genes=function(genes){
  ret=unlist(strsplit(genes,",",fixed=T))
  ret=sapply(ret,function(x){
    x=gsub(" ","",x)
  })
  return(ret)
}

r=function() runApp(launch.browser=F,port=8884)
