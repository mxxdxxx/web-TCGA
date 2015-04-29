get.Data.Methylation.Heatmap = function(input){
  
  my.Genes = input$methylation.Genes
  my.Entity = input$methylation.Entity
  
  if(input$methylation.P.or.Diff=="P-Value" & input$methylation.Paired.Only){
    my.Fun = function(x,y) return(wilcox.test(x,y,paired=T)$p.value)
  }else{
    my.Fun = function(x,y) return(wilcox.test(x,y,paired=F)$p.value)
  }
  if(input$methylation.P.or.Diff=="Difference"){
    my.Fun = function(x,y) return(mean(x, na.rm=T) - mean(y, na.rm=T))
  }
  
  methylation.Table = methylation.List[[my.Entity]]
  
  normal.Pats = get.Normals(unique(methylation.Table$patient))
  tumor.Pats = get.Tumors(unique(methylation.Table$patient))
  normal.Pats = normal.Pats[which(getTSS(normal.Pats) == my.Entity)]
  tumor.Pats = tumor.Pats[which(getTSS(tumor.Pats) == my.Entity)]
  pairs = get.Pairs(c(normal.Pats, tumor.Pats))
  
  if(length(pairs) == 0 & input$methylation.Paired.Only){
    return("paired-error")
  }
  
  if(input$methylation.Paired.Only){
    to.Join = pos.2.gene[J(my.Genes)][, list(IlmnID, UCSC_RefGene_Name, UCSC_RefGene_Group)]
    setkey(to.Join, "IlmnID")
    meth.Subset = methylation.Table[to.Join]
    meth.Subset = meth.Subset[which(meth.Subset$patient %in% unlist(pairs)), ]
    meth.Subset = cbind(meth.Subset[is.Tumor(meth.Subset$patient)], meth.Subset[!is.Tumor(meth.Subset$patient), beta])
    
    meth.Subset[, val:=lapply(.SD, function(x=beta,y=V2) my.Fun(x,y)),
                by=list(patient,
                        UCSC_RefGene_Name,
                        UCSC_RefGene_Group),
                .SDcols=c("beta","V2") ]
    meth.Subset[, rn:=NULL]
    meth.Subset[, V2:=NULL]
    meth.Subset[, beta:=NULL]
    meth.Subset$patient = substr(meth.Subset$patient, 1, 12)
    meth.Subset = meth.Subset[!duplicated(meth.Subset)]
    to.Plot = meth.Subset
  }else{
    
    if(length(normal.Pats) == 0){
      return("paired-error")
    }else{
      
      to.Join = pos.2.gene[J(my.Genes)][, list(IlmnID, UCSC_RefGene_Name, UCSC_RefGene_Group)]
      setkey(to.Join, "IlmnID")
      meth.Subset = methylation.Table[to.Join]
      tumor.Subset = meth.Subset[which(meth.Subset$patient %in% unlist(tumor.Pats)), ]
      
      # Generate normal list by, region
      normal.Subset = meth.Subset[which(meth.Subset$patient %in% unlist(normal.Pats))]
      normal.List = lapply(unique(normal.Subset$UCSC_RefGene_Name), function(gene){
        to.Name = lapply(unique(normal.Subset$UCSC_RefGene_Group), function(region){
          normal.Subset[normal.Subset$UCSC_RefGene_Name==gene & normal.Subset$UCSC_RefGene_Group==region,][, beta]
          })
        names(to.Name) = unique(normal.Subset$UCSC_RefGene_Group)
        to.Name
        })
      names(normal.List) = unique(normal.Subset$UCSC_RefGene_Name)
            
      tumor.Subset[, val:=lapply(.SD, function(x=beta,y=UCSC_RefGene_Name,z=UCSC_RefGene_Group){
        y=unique(y)
        z=unique(z)
        if(!is.numeric(x) | !y %in% names(normal.List) | !z %in% names(normal.List[[y]]) ) return(NaN)
        my.Fun(x, normal.List[[y]][[z]])
        }),
        by=list(patient,
                UCSC_RefGene_Name,
                UCSC_RefGene_Group),
        .SDcols=c("beta","UCSC_RefGene_Name", "UCSC_RefGene_Group")]
      
      tumor.Subset$patient = substr(tumor.Subset$patient, 1, 12)
      tumor.Subset[, rn:=NULL]
      tumor.Subset[, beta:=NULL]
      tumor.Subset = tumor.Subset[!duplicated(tumor.Subset)]
      
      to.Plot = tumor.Subset
    }
  }
  
  setnames(to.Plot, "UCSC_RefGene_Name", "Gene")
  setnames(to.Plot, "UCSC_RefGene_Group", "Part")
  
  # Correct p-Values
  if(input$methylation.P.or.Diff=="P-Value"){
    setattr(to.Plot, "val", p.adjust(to.Plot$val, "hochberg"))
    to.Plot[, is.Differential := to.Plot$val < 0.05]
    to.Join = to.Plot[,(sum(is.Differential) /length(is.Differential))*100 ,by=list(Gene,Part)]
    to.Join.2 = to.Plot[,sum(is.Differential) ,by=list(Gene,Part)]
    setkey(to.Join, "Gene", "Part")
    setkey(to.Join.2, "Gene", "Part")
    setkey(to.Plot, "Gene", "Part")
    to.Plot = to.Plot[to.Join]
    to.Plot = to.Plot[to.Join.2]
    setnames(to.Plot, "V1", "Percentage")
    setnames(to.Plot, "V1.1", "p.Amount")
    
  }else{
    to.Plot[, is.Differential := (to.Plot$val > input$methylation.Threshold | to.Plot$val < (input$methylation.Threshold*-1))]
    to.Join = to.Plot[,(sum(is.Differential) /length(is.Differential))*100 ,by=list(Gene,Part)]
    to.Join.2 = to.Plot[,sum(is.Differential) ,by=list(Gene,Part)]
    setkey(to.Join, "Gene", "Part")
    setkey(to.Join.2, "Gene", "Part")
    setkey(to.Plot, "Gene", "Part")
    to.Plot = to.Plot[to.Join]
    to.Plot = to.Plot[to.Join.2]
    setnames(to.Plot, "V1", "d.Percentage")
    setnames(to.Plot, "V1.1", "d.Amount")
  }
  return(as.data.frame(to.Plot))
}


plot.Methylation.Heatmap = function(input){
  
  to.Plot = get.Data.Methylation.Heatmap(input)
  
  if(to.Plot == "paired-error")
    return(plot.new()+mtext("This entity does not have paired samples", col="red"))

  d = ggplot(to.Plot, aes(x = Gene, y = Part))
  
  if(input$methylation.P.or.Diff == "P-Value" & input$methylation.Use.Percentage){
    d = d + stat_sum(aes(size = Percentage))
  }
  if(input$methylation.P.or.Diff == "P-Value" & !input$methylation.Use.Percentage){
    d = d + stat_sum(aes(size = p.Amount))
  }
  if(input$methylation.P.or.Diff == "P-Difference" & input$methylation.Use.Percentage){
    d = d + stat_sum(aes(size = d.Percentage))
  }
  if(input$methylation.P.or.Diff == "Difference" & !input$methylation.Use.Percentage){
    d = d + stat_sum(aes(size = d.Amount))
  }
  return(d)
}

plot.Methylation.Distribution = function(input, i){
  
  if(input$methylation.P.or.Diff!="Difference")
    return(plot.new()+mtext("Please choose \"Difference for Calculation\"", col="red"))
  
  my.Gene = input$methylation.Genes[i]
  to.Plot = get.Data.Methylation.Heatmap(input)
  to.Plot = to.Plot[to.Plot$Gene == my.Gene, ]
  
  if(to.Plot == "paired-error")
    return(plot.new()+mtext("This entity does not have paired samples", col="red"))
  
  d = ggplot(to.Plot, aes(x=val))
  d = d + geom_histogram() +
    ggtitle(paste("Distribution of ", my.Gene, "s Differential Methylation", sep=""))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    xlab("Difference of Tumor Methylation")
  return(d)
}

plot.Methylation.Waterfall = function(input, i){
  
  if(input$methylation.P.or.Diff!="Difference")
    return(plot.new()+mtext("Please choose \"Difference for Calculation\"", col="red"))
  my.Gene = input$methylation.Genes[i]
  to.Plot = get.Data.Methylation.Heatmap(input)
  to.Plot = to.Plot[to.Plot$Gene == my.Gene, ]
  to.Plot = to.Plot[which(to.Plot$is.Differential)]
  if(dim(to.Plot)[1]==0) return(plot.new()+mtext("This Gene has no differential methylated Patients", col="red"))
  
  if(to.Plot == "paired-error")
    return(plot.new() + mtext("This entity does not have paired samples", col="red"))
  
  to.Plot = to.Plot[order(to.Plot$val, decreasing=T), ]  
  to.Plot = cbind(to.Plot, id=1:nrow(to.Plot))
  to.Plot = cbind(to.Plot, start=replicate(nrow(to.Plot), 0))
  to.Plot = cbind(to.Plot,x.Axis.Id = paste(to.Plot$patient,
                                            to.Plot$Gene,
                                            to.Plot$Part,
                                            to.Plot$Difference,
                                            sep="_"))
  
  d = ggplot(to.Plot, aes(x.Axis.Id, fill = Part, label=to.Plot$Tumor),
             environment = environment()) 
  d = d +
    geom_rect(
      aes(x = x.Axis.Id,
          xmin = id ,
          xmax = id + 0.75,
          ymin = val,
          ymax = start)) +
    theme(axis.text.x=element_text(angle=90, hjust=1, size=5)) +
    ggtitle(paste("Methylation of", input$my.Gene )) +
    scale_x_discrete("Patients", labels=as.character(to.Plot$patient))
  return(d)
}

plot.Methylation.Table = function(x){
  to.Plot = get.Data.Methylation.Heatmap(input)
  return(to.Plot)
}