getTSS=function(barcodes){
  pats.tss=sapply(barcodes, function(b){
    if(b != "rn"){
      tss=substr(b,6,7)
      TSS[TSS$TSS.Code==tss,]$Study.Name
    } else {
      "Index"
    }
  })
  return(pats.tss)
}

getCancerTypes=function(pat_ids){
  pat_ids=setdiff(pat_ids,"rn")
  res=getTSS(pat_ids)
  res=unique(unlist(res))
  return(res)
}

get.Pairs = function(barcodes){
  
  patients = sapply(barcodes, function(barcode){
    if(grepl("tcga", barcode, ignore.case=T)){
      substr(barcode, 1, 12)
    }else{
      NA
    }
  })
  pairs = table(patients)
  pairs = names(pairs[which(pairs >= 2)])
  paired.Barcodes = lapply(pairs, function(patient) barcodes[grep(patient, barcodes)])
  
  paired.Barcodes = lapply(paired.Barcodes, function(pair){
    if(is.Tumor(pair[1]) & is.Normal(pair[2]))
      pair
    if(is.Tumor(pair[2]) & is.Normal(pair[1]))
      pair
    })
  paired.Barcodes = paired.Barcodes[!sapply(paired.Barcodes, is.null)]
  
  return(paired.Barcodes)
}

get.Normals = function(barcodes){
  ret = sapply(barcodes, function(barcode){
    if(substr(barcode, 14,14) == 1)
      barcode
    else
      NULL
    })
  names(ret)[!sapply(ret, is.null)]
}

get.Tumors = function(barcodes){
  ret = sapply(barcodes, function(barcode){
    if(substr(barcode, 14,14) == 0)
      barcode
    else
      NULL
  })
  names(ret)[!sapply(ret, is.null)]
}

is.Tumor = function(x) return(substr(x, 14,14) == 0)

is.Normal = function(x) return(substr(x, 14,14) == 1)

has.Pair = function(sample, barcodes){
  res = table(substr(sample, 1, 12) == substr(barcodes, 1, 12))
  if(length(res)>1 & res["TRUE"]>1)
    return(T)
  else
    return(F)
}

get.Pair = function(sample, barcodes){
  pair = which(substr(sample, 1, 12) == substr(barcodes, 1, 12))
  return(barcodes[pair[which(barcodes[pair] != sample)]])
}





