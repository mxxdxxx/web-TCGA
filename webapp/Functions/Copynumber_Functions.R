generate.Copynumber.Matrix = function(genes, entities){
  ret = list()
  i=0
  for(g in genes){
    for(e in entities){
      i = i+1
      pats = which(getTSS(colnames(copynumber.Table)) == e)
      mat = unlist(subset(copynumber.Table[J(g)], T, pats))
      res = table(c(mat, -2, -1 , 0, 1, 2))
      res = res -1
      res = c(sum(res["-2"], res["2"]), sum(res["-1"], res["1"]),
              sum(res["-2"], res["2"], res["-1"], res["1"]))
      res = t(do.call(cbind, replicate(length(mat), res, simplify=F)))
      res = cbind(mat, res, round((res/length(mat))*100, 1))
      colnames(res)[2:ncol(res)] = c("High", "Low", "Both",
                                     "High.P", "Low.P", "Both.P")
      res = cbind(res, Gene = replicate(nrow(res),g),
                  Entity = replicate(nrow(res),e))
      ret[[i]] = res
    }
  }
  ret = as.data.frame(do.call(rbind, ret))
  for(cc in c("High", "Low", "Both","High.P", "Low.P", "Both.P")){
    ret[, cc] = as.numeric(as.character(ret[, cc]))
  }
  return(ret)
}


plot.Global.Copynumber = function(input){
  
  genes = input$copynumber.Genes
  entities = input$copynumber.Entities
  
  to.Plot = generate.Copynumber.Matrix(genes, entities)
  
  y.axis = element_text(color = "black")
  x.axis = element_text(face = "italic",color = "black")

  if(input$copynumber.Hi.Low == "Both" & input$copynumber.Use.Percentage){
    d = ggplot(to.Plot, aes(x=Gene, y=Entity)) +
      stat_sum(aes(size = Both.P )) +
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      ggtitle("Global Copy Number Variation Status") +
      ylab("Tumor Entity") +
      xlab("Genes") +
      labs(size = "Gain/Loss\nRate (%)") +
      theme(axis.text.x = x.axis, axis.text.y = y.axis)
  }
  if(input$copynumber.Hi.Low == "Low Level Only" & input$copynumber.Use.Percentage){
    d = ggplot(to.Plot, aes(x=Gene, y=Entity)) +
      stat_sum(aes(size=Low.P)) +
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      ggtitle("Global Copy Number Variation Status") +
      ylab("Tumor Entity") +
      xlab("Genes") +
      labs(size = "Gain/Loss\nRate (%)") +
      theme(axis.text.x = x.axis, axis.text.y = y.axis)
  }
  if(input$copynumber.Hi.Low == "High Level Only" & input$copynumber.Use.Percentage){
    d = ggplot(to.Plot, aes(x=Gene, y=Entity)) +
      stat_sum(aes(size=High.P)) +
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      ggtitle("Global Copy Number Variation Status") +
      ylab("Tumor Entity") +
      xlab("Genes") +
      labs(size = "Gain/Loss\nRate (%)") +
      theme(axis.text.x = x.axis, axis.text.y = y.axis)
  }
  if(input$copynumber.Hi.Low == "Both" & !input$copynumber.Use.Percentage){
    d = ggplot(to.Plot, aes(x=Gene, y=Entity)) +
      stat_sum(aes(size=Both)) +
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      ggtitle("Global Copy Number Variation Status") +
      ylab("Tumor Entity") +
      xlab("Genes") +
      labs(size = "Gain/Loss\nRate (%)") +
      theme(axis.text.x = x.axis, axis.text.y = y.axis)
  }
  if(input$copynumber.Hi.Low == "Low Level Only" & !input$copynumber.Use.Percentage){
    d = ggplot(to.Plot, aes(x=Gene, y=Entity)) +
      stat_sum(aes(size=Low)) +
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      ggtitle("Global Copy Number Variation Status") +
      ylab("Tumor Entity") +
      xlab("Genes") +
      labs(size = "Gain/Loss\nRate (%)") +
      theme(axis.text.x = x.axis, axis.text.y = y.axis)
  }
  if(input$copynumber.Hi.Low == "High Level Only" & !input$copynumber.Use.Percentage){
    d = ggplot(to.Plot, aes(x=Gene, y=Entity)) +
      stat_sum(aes(size=High)) +
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      ggtitle("Global Copy Number Variation Status") +
      ylab("Tumor Entity") +
      xlab("Genes") +
      labs(size = "Gain/Loss\nRate (%)") +
      theme(axis.text.x = x.axis, axis.text.y = y.axis)
  }
  
  return(d)
}


plot.Copynumber.Distribution = function(input, i){
  
  genes = input$copynumber.Genes
  entity = input$copynumber.Entities[i]
  
  to.Plot = generate.Copynumber.Matrix(genes, entity)
  to.Plot = to.Plot[to.Plot$Entity == entity, ]
  
  cols = c("-2" = "#ca0020", "-1" = "#f4a582", "0" = "#f7f7f7", "1" = "#92c5de",
           "2" = "#0571b0")
  
  d = ggplot(to.Plot, aes(Gene, fill=mat))
  d = d + geom_bar() +
    theme(axis.text.x=element_text(face="italic")) +
    guides(fill=guide_legend(title="CNV Level")) +
    ylab("Patient Count") +
    xlab("Genes") +
    ggtitle(paste("Copy Number Variation Status of", entity)) +
    scale_fill_manual(values=cols, breaks=c("-2","-1","0","1","2"))
  return(d)
}
