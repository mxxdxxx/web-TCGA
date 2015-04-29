
# Get Table Representation ------------------------------------------------

variant.Calculate.Table=function(input){

  my.Genes = input$variant.Genes
  my.Entities = input$variant.Entities
  
  if(input$variant.Exclude.Silent){
    local.Variant.Table.silent = variant.Table
    local.Variant.Table = subset(variant.Table, Variant_Classification!="Silent")
  } else {
    local.Variant.Table = variant.Table
  }

  if(length(my.Genes)==0 | length(my.Entities)==0)
    return(plot.new()+mtext("No Genes and/or Cancer Types entered", col="red"))
  
  local.Variant.Table = subset(local.Variant.Table,
                               Hugo_Symbol %in% my.Genes &
                                 Entity %in% my.Entities)
  
  if(input$variant.Multi.Mut.Per.Gene){
    snp.Pat = matrix(nrow=length(unique(local.Variant.Table$Hugo_Symbol)),
                     ncol=length(unique(local.Variant.Table$Entity)),
                     dimnames=list(unique(local.Variant.Table$Hugo_Symbol),
                                   unique(local.Variant.Table$Entity)))
    for(g in rownames(snp.Pat)){
      for(t in colnames(snp.Pat)){
        mt = subset(local.Variant.Table[J(g)], Entity == t )
        all.Pats = unlist(subset(variant.Table, Entity==t ,
                          "Tumor_Sample_Barcode" ))
        snp.Pat[g,t] = (nrow(mt) / length(unique(all.Pats)))*100
      }
    }
    local.Variant.Table = cbind(local.Variant.Table,rep(0,nrow(local.Variant.Table)))
    setnames(local.Variant.Table, ncol(local.Variant.Table), "mut_rate")
    for(i in 1:nrow(local.Variant.Table)){
      set(local.Variant.Table, i, ncol(local.Variant.Table),
          snp.Pat[local.Variant.Table$Hugo_Symbol[i],
                  as.character(local.Variant.Table$Entity[i])])
    }
    
    return(local.Variant.Table)
    
  } else {
    snp.Pat = matrix(nrow=length(unique(local.Variant.Table$Hugo_Symbol)),
                     ncol=length(unique(local.Variant.Table$Entity)),
                     dimnames=list(unique(local.Variant.Table$Hugo_Symbol),
                                   unique(local.Variant.Table$Entity)))
    for(g in rownames(snp.Pat)){
      for(t in colnames(snp.Pat)){
        mt = subset(local.Variant.Table[J(g)], Entity == t )
        all.Pats = unlist(subset(variant.Table, Entity==t ,
                                 "Tumor_Sample_Barcode" ))
        snp.Pat[g,t] = (length(unique(mt$Tumor_Sample_Barcode)) /
                          length(unique(all.Pats)))*100
        
      }
    }
    local.Variant.Table = cbind(local.Variant.Table,
                                rep(0,nrow(local.Variant.Table)))
    setnames(local.Variant.Table, ncol(local.Variant.Table), "mut_rate")
    for(i in nrow(local.Variant.Table)){
      set(local.Variant.Table, i, ncol(local.Variant.Table),
          snp.Pat[local.Variant.Table$Hugo_Symbol[i],
                  as.character(local.Variant.Table$Entity[i])])
    }
    return(local.Variant.Table)
  }
}


# Pie Charts --------------------------------------------------------------

plot.Variant.PieCharts=function(input,i){
  
  my.Genes=input$variant.Genes
  my.Entity=input$variant.Entities[i]
  
  if(input$variant.Exclude.Silent){
    local.Variant.Table.silent = variant.Table
    local.Variant.Table = subset(variant.Table, Variant_Classification!="Silent")
  } else {
    local.Variant.Table = variant.Table
  }
  
  if(length(my.Genes)==0 | length(my.Entity)==0)
    return(plot.new()+mtext("No Genes and/or Cancer Types entered", col="red"))
  
  all.Pats = unlist(subset(local.Variant.Table, Entity==my.Entity,
                           "Tumor_Sample_Barcode"))
  
  if(input$variant.Multi.Mut.Per.Gene){
    
    n.Pats = length(unique(all.Pats))
    local.Variant.Table = subset(local.Variant.Table,
                                 Hugo_Symbol %in% my.Genes &
                                   Entity %in% my.Entity )
    if(input$mut_type_split){
      to.Pie = table(as.character(local.Variant.Table$Combined_Class))
    }else{
      to.Pie = table(as.character(local.Variant.Table$Hugo_Symbol))
    }
  }else{
    n.Pats = length(unique(all.Pats))
    local.Variant.Table = subset(local.Variant.Table,
                                 Hugo_Symbol %in% my.Genes &
                                   Entity %in% my.Entity )
    marc.Duplicates = subset(local.Variant.Table, T,
                             c(Hugo_Symbol,Tumor_Sample_Barcode))
    local.Variant.Table = subset(local.Variant.Table,
                                 !duplicated(marc.Duplicates))
    if(input$mut_type_split){
      to.Pie = table(as.character(local.Variant.Table$Combined_Class))
    }else{
      to.Pie = table(as.character(local.Variant.Table$Hugo_Symbol))
    }
  }
  
  to.Pie=to.Pie[which( (to.Pie/n.Pats)*100 > input$threshold)]
  
  return(
    pie(to.Pie,
        main=paste("Mutations of ",my.Entity,sep=""),
        cex=.65,
        clockwise=T,
        init.angle=90,
        labels=paste(round(to.Pie/n.Pats,digits=3)*100,"% ",gsub("[[:punct:]]"," ", names(to.Pie)),sep="")
    )
  )
}


# Heatmap -----------------------------------------------------------------

plot.Variant.Heatmap = function(input){
    
  to.Plot  = variant.Calculate.Table(input)
    
  y.axis = element_text(color = "black")
  x.axis = element_text(face = "italic",color = "black")
  ret = ggplot(to.Plot, aes(x=Hugo_Symbol, y=Entity)) +
    stat_sum(aes(size=mut_rate)) +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    ggtitle("Mutational Landscape") +
    ylab("Tumor Entity") +
    xlab("Genes") +
    labs(size = "Mutation\nRate (%)") +
    theme(axis.text.x = x.axis,axis.text.y = y.axis)
  return(ret)
  
}

