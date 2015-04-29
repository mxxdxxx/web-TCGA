folchange = function(a,b){
  a/b
}

# Subset z Score Matrix ---------------------------------------------------

generate.Plot.Matrix = function(genes, entities, threshold){
  
  if(missing(threshold)) threshold = 2
  
  tss = getTSS(colnames(expression.Table))
  
  entity.2.pat = lapply(entities, function(x){
    pats = names(tss)[which(tss == x)]
    tumors = get.Tumors(pats)
  })
  names(entity.2.pat) = entities
  
  expression.By.Entity = list()
  for(pats in names(entity.2.pat)){
    pat.Ids = entity.2.pat[[pats]]
    expression.By.Gene = list()
    for(gene in genes){
      tmp.Dt = unlist(subset(expression.Table[J(gene)], T, pat.Ids))
      res = lapply(1:length(tmp.Dt), function(x){
        exp = tmp.Dt[x]
        rest = tmp.Dt[c(1:(x-1),(x+1):10)]
        rest.Mean = mean(rest)
        z.Score = (exp - rest.Mean) / sd(rest)
        fold.Change = exp / rest.Mean
        c(z.Score, fold.Change)
      })
      res = do.call(rbind, res)
      res = as.data.frame(res)
      rownames(res) = pat.Ids
      res = cbind(res, replicate(nrow(res), gene), replicate(nrow(res), pats))
      colnames(res) = c("z.Score", "fold.Change", "Gene", "Entity")
      z.Exp = length(which(res$z.Score > threshold | res$z.Score < (threshold*-1)))
      fold.Exp = length(which(res$fold.Change > threshold | res$fold.Change < (threshold*-1)))
      
      res = cbind(res,
                  z.Amount = z.Exp,
                  fold.Amount = fold.Exp,
                  z.Percentage = (z.Exp / nrow(res))*100,
                  fold.Percentage = (fold.Exp / nrow(res))*100)
      
      expression.By.Gene[[gene]] = res
    }
    expression.By.Entity[[pats]] = do.call(rbind, expression.By.Gene)
    expression.By.Entity[[pats]] = cbind(expression.By.Entity[[pats]],
                                         Patient = pat.Ids)
  }
  to.Plot = do.call(rbind, expression.By.Entity)
  
  return(to.Plot)
}


plot.Global.Expresion = function(input){
  
  to.Plot = generate.Plot.Matrix(input$expression.Genes,
                                 input$expression.Entities,
                                 threshold = input$expression.Threshold)
  
  
  if(input$expression.Z.or.Foldchange == "Z-Score" & input$expression.Use.Percentage){
    d = ggplot(to.Plot, aes(x = Entity, y = Gene))
    d = d +
      stat_sum(aes(colour = Entity, size = z.Percentage)) +
      ggtitle("Number of Over or Under Expressed Genes") +
      theme(axis.text.x=element_text(angle=90, hjust=1),
            axis.text.y=element_text(hjust=1, face="italic"))+
      labs(size="Percentage")
  }
  if(input$expression.Z.or.Foldchange == "Z-Score" & !input$expression.Use.Percentage){
    d = ggplot(to.Plot, aes(x = Entity, y = Gene))
    d = d +
      stat_sum(aes(colour = Entity, size = z.Amount)) +
      ggtitle("Number of Over or Under Expressed Genes") +
      theme(axis.text.x=element_text(angle=90, hjust=1),
            axis.text.y=element_text(hjust=1, face="italic"))+
      labs(size="#Patients")
  }
  if(input$expression.Z.or.Foldchange == "Foldchange" & input$expression.Use.Percentage){
    d = ggplot(to.Plot, aes(x = Entity, y = Gene))
    d = d +
      stat_sum(aes(colour = Entity, size = fold.Percentage)) +
      ggtitle("Number of Over or Under Expressed Genes") +
      theme(axis.text.x=element_text(angle=90, hjust=1),
            axis.text.y=element_text(hjust=1, face="italic"))+
      labs(size="Percentage")
  }
  if(input$expression.Z.or.Foldchange == "Foldchange" & !input$expression.Use.Percentage){
    d = ggplot(to.Plot, aes(x = Entity, y = Gene))
    d = d +
      stat_sum(aes(colour = Entity, size = fold.Amount)) +
      ggtitle("Number of Over or Under Expressed Genes") +
      theme(axis.text.x=element_text(angle=90, hjust=1),
            axis.text.y=element_text(hjust=1, face="italic"))+
      labs(size="#Patients")
  }
  return(d)
}


plot.Expression.Waterfall.by.Gene = function(input, i){
  
  to.Plot = generate.Plot.Matrix(input$expression.Genes[i],
                                 input$expression.Entities,
                                 threshold = input$expression.Threshold)
  
  if(input$expression.Z.or.Foldchange == "Z-Score"){
    to.Keep = which(to.Plot$z.Score > input$expression.Threshold |
                      to.Plot$z.Score < (input$expression.Threshold*-1))
    to.Plot = to.Plot[to.Keep, ]
    to.Plot = to.Plot[order(to.Plot$z.Score, decreasing=T), ]
    
    to.Plot = cbind(to.Plot, id=1:nrow(to.Plot))
    to.Plot = cbind(to.Plot, start=replicate(nrow(to.Plot), 0))
    to.Plot = cbind(to.Plot,x.Axis.Id = paste(to.Plot$Patient, to.Plot$Gene,
                                              to.Plot$z.Score, sep="_"))
    
    d = ggplot(to.Plot, aes(x.Axis.Id, fill = Entity, label=to.Plot$Patient),
               environment = environment())
    d = d +
      geom_rect(
        aes(x = x.Axis.Id,
            xmin = id ,
            xmax = id + 0.75,
            ymin = z.Score,
            ymax = start)) +
      theme(axis.text.x=element_text(angle=90, hjust=1, size=5)) +
      ggtitle(paste("Expression of ", input$expression.Genes[i], ", n = ",
                    length(unique(to.Plot$Patient)), sep="")) +
      scale_x_discrete("Patients", labels=as.character(to.Plot$Patient))
  } else {
    to.Keep = which(to.Plot$fold.Change > input$expression.Threshold |
                      to.Plot$fold.Change < (input$expression.Threshold*-1))
    to.Plot = to.Plot[to.Keep, ]
    to.Plot = to.Plot[order(to.Plot$fold.Change, decreasing=T), ]
    
    to.Plot = cbind(to.Plot, id=1:nrow(to.Plot))
    to.Plot = cbind(to.Plot, start=replicate(nrow(to.Plot), 0))
    to.Plot = cbind(to.Plot,x.Axis.Id = paste(to.Plot$Patient, to.Plot$Gene,
                                              to.Plot$fold.Change, sep="_"))
    
    d = ggplot(to.Plot, aes(x.Axis.Id, fill = Entity, label=to.Plot$Patient),
               environment = environment())
    d = d +
      geom_rect(
        aes(x = x.Axis.Id,
            xmin = id ,
            xmax = id + 0.75,
            ymin = fold.Change,
            ymax = start)) +
      theme(axis.text.x=element_text(angle=90, hjust=1, size=5)) +
      ggtitle(paste("Expression of", input$expression.Genes[i] )) +
      scale_x_discrete("Patients", labels=as.character(to.Plot$Patient))
  }
  return(d)
}


plot.Expression.Waterfall.by.Entity = function(input, i){
  
  to.Plot = generate.Plot.Matrix(input$expression.Genes,
                                 input$expression.Entities[i],
                                 threshold = input$expression.Threshold)
  
  if(input$expression.Z.or.Foldchange == "Z-Score"){
    to.Keep = which(to.Plot$z.Score > input$expression.Threshold |
                      to.Plot$z.Score < (input$expression.Threshold*-1))
    to.Plot = to.Plot[to.Keep, ]
    to.Plot = to.Plot[order(to.Plot$z.Score, decreasing=T), ]
    
    to.Plot = cbind(to.Plot, id=1:nrow(to.Plot))
    to.Plot = cbind(to.Plot, start=replicate(nrow(to.Plot), 0))
    to.Plot = cbind(to.Plot,x.Axis.Id = paste(to.Plot$Patient,
                                              to.Plot$Gene,
                                              to.Plot$Entity
                                              , sep="_"))
    
    d = ggplot(to.Plot, aes(x.Axis.Id, fill = Gene, label=to.Plot$Patient),
               environment = environment())
    d = d +
      geom_rect(
        aes(x = x.Axis.Id,
            xmin = id ,
            xmax = id + 0.75,
            ymin = z.Score,
            ymax = start)) +
      theme(axis.text.x=element_text(angle=90, hjust=1, size=5)) +
      ggtitle(paste("Expression of", input$expression.Entities[i] )) +
      scale_x_discrete("Patients", labels=as.character(to.Plot$Patient))
  } else {
    to.Keep = which(to.Plot$fold.Change > input$expression.Threshold |
                      to.Plot$fold.Change < (input$expression.Threshold*-1))
    to.Plot = to.Plot[to.Keep, ]
    to.Plot = to.Plot[order(to.Plot$fold.Change, decreasing=T), ]
    
    to.Plot = cbind(to.Plot, id=1:nrow(to.Plot))
    to.Plot = cbind(to.Plot, start=replicate(nrow(to.Plot), 0))
    to.Plot = cbind(to.Plot,x.Axis.Id = paste(to.Plot$Patient, to.Plot$Gene,
                                              to.Plot$fold.Change, sep="_"))
    
    d = ggplot(to.Plot, aes(x.Axis.Id, fill = Entity, label=to.Plot$Patient),
               environment = environment())
    d = d +
      geom_rect(
        aes(x = x.Axis.Id,
            xmin = id ,
            xmax = id + 0.75,
            ymin = fold.Change,
            ymax = start)) +
      theme(axis.text.x=element_text(angle=90, hjust=1, size=5)) +
      ggtitle(paste("Expression of", input$expression.Genes[i] )) +
      scale_x_discrete("Patients", labels=as.character(to.Plot$Patient))
  }
  return(d)
}


plot.Expression.By.Entity = function(input, i){
  
  to.Plot = generate.Plot.Matrix(input$expression.Genes,
                                     input$expression.Entities[i])
  
  if(input$expression.Z.or.Foldchange == "Z-Score"){
    d = ggplot(to.Plot, aes(x=Gene, y=z.Score))
    d = d +
      geom_boxplot(outlier.colour = "red", outlier.size = 3) +
      ggtitle(paste("Relative Expression of", input$expression.Entities[i] )) +
      ylab("Z-Score") +
      theme(axis.text.x=element_text(angle=90, hjust=1, face="italic"))
  }else{
    d = ggplot(to.Plot, aes(x=Gene, y=fold.Change))
    d = d +
      geom_boxplot(outlier.colour = "red", outlier.size = 3) +
      ggtitle(paste("Relative Expression of", input$expression.Entities[i] )) +
      ylab("Fold Change") +
      theme(axis.text.x=element_text(angle=90, hjust=1, face="italic"))  
  }
  
  return(d)
}


plot.Expression.By.Gene = function(input, i){  
  
  to.Plot = generate.Plot.Matrix(input$expression.Genes[i],
                                 input$expression.Entities)
  
  if(input$expression.Z.or.Foldchange == "Z-Score"){
    d = ggplot(to.Plot, aes(x=Entity, y=z.Score))
    d = d +
      geom_boxplot(outlier.colour = "red", outlier.size = 3) +
      ggtitle(paste("Relative Expression of", input$expression.Genes[i] )) +
      ylab("Z-Score") +
      theme(axis.text.x=element_text(angle=90, hjust=1, face="italic"))
  }else{
    d = ggplot(to.Plot, aes(x=Entity, y=fold.Change))
    d = d +
      geom_boxplot(outlier.colour = "red", outlier.size = 3) +
      ggtitle(paste("Relative Expression of", input$expression.Genes[i] )) +
      ylab("Fold Change") +
      theme(axis.text.x=element_text(angle=90, hjust=1, face="italic"))
  }
  return(d)
}
