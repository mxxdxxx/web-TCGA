max_plots=50
hw=750

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output, session) {
  
  
  # Variant Section ---------------------------------------------------------  
  
  #Header
  output$caption=renderText({
    "Heatmap of your selected genes and cancer types"
  })
  
  #variant.Table
  output$variant.Table = renderDataTable({
    drops=c("Combined","Combined_Class","End")
    to.tab=as.data.frame(variant.Calculate.Table(input))
    to.tab=to.tab[,which(!colnames(to.tab) %in% drops)]
    cn=gsub("[[:punct:]]",replacement=" ",x=colnames(to.tab))
    cn=gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", cn, perl=TRUE)
    colnames(to.tab)=cn
    to.tab
  })
  
  
  # Variant Reactive plotting Section, create the heatmap -------------------  
  output$plot.Variant.Heatmap = renderPlot({
    withProgress(session, {
      setProgress(message = "Creating Plots",
                  detail = "This may take a few seconds...")
      print(plot.Variant.Heatmap(input))
    })
  })
  
  #Create PIE-charts
  output$plot.Variant.PieCharts = renderUI({
    plot_output_list = lapply(1:length(input$variant.Entities), function(i) {
      plotname = paste("variant.PieChart", i, sep="")
      plotOutput(plotname, height = hw, width = hw)
    })
    do.call(tagList, plot_output_list)
  })
  
  for (i in 1:max_plots) {
    local({
      my_i = i
      plotname = paste("variant.PieChart", my_i, sep="")
      output[[plotname]] = renderPlot({
        withProgress(session, {
          setProgress(message = "Creating Plots",
                      detail = "This may take a few seconds...")
          plot.Variant.PieCharts(input,my_i)
        })
      })
    })
  }
  
  
  # Variant Section, Downloads ----------------------------------------------  
  output$download.Plot.Variant.Heatmap = downloadHandler(
    filename = "Variant_Heatmap.pdf",
    content = function(file) {
      pdf(file)
      print(plot.Variant.Heatmap(input))
      dev.off()
    }
  )
  
  #handling the PIE Charts
  output$download.Plot.Variant.PieCharts = downloadHandler(
    filename = "Variant_Pie-Charts.pdf",
    content = function(file) {
      pdf(file)
      plot_output_list = lapply(1:length(input$cancer_types), function(i) {
        plotname = paste("plot", i, sep="")
        print(plot.Variant.PieCharts(input,i))
      })
      dev.off()
    }
  )
  
  #serve data for variant.Table download
  output$download.Variant.Table=downloadHandler(
    filename = "Variants.csv",
    content = function(file) {
      write.csv(variant.Calculate.Table(input),file,sep=",")
    }
  )
  
  
  # Methylation Section, Heatmap -----------------------------------------
  
  output$plot.Methylation.Heatmap=renderPlot({
    withProgress(session, {
      setProgress(message = "Creating Plots",
                  detail = "This may take a few seconds...")
      print(plot.Methylation.Heatmap(input))
    })
  })
  
    
  output$plot.Methylation.Distribution = renderUI({
    number_of_plots=length(input$methylation.Genes)
    plot_output_list = lapply(1:number_of_plots, function(i) {
      plotname = paste("plot.Methylation.Distribution", i, sep="")
      plotOutput(plotname, height = hw, width = hw)
    })
    do.call(tagList, plot_output_list)
  })
  
  for (i in 1:max_plots) {
    local({
      my_i = i
      plotname = paste("plot.Methylation.Distribution", my_i, sep="")
      output[[plotname]] = renderPlot({
        withProgress(session, {
          setProgress(message = "Creating Plots",
                      detail = "This may take a few seconds...")
          print(plot.Methylation.Distribution(input,my_i))
        })
      })
    })
  }
  
  output$plot.Methylation.Waterfall = renderUI({
    number_of_plots=length(input$methylation.Genes)
    plot_output_list = lapply(1:number_of_plots, function(i) {
      plotname = paste("plot.Methylation.Waterfall", i, sep="")
      plotOutput(plotname, height = hw, width = hw)
    })
    do.call(tagList, plot_output_list)
  })
  
  for (i in 1:max_plots) {
    local({
      my_i = i
      plotname = paste("plot.Methylation.Waterfall", my_i, sep="")
      output[[plotname]] = renderPlot({
        withProgress(session, {
          setProgress(message = "Creating Plots",
                      detail = "This may take a few seconds...")
          print(plot.Methylation.Waterfall(input,my_i))
        })
      })
    })
  }
  
  output$methylation.Table = renderDataTable({
    get.Data.Methylation.Heatmap(input)
  })
    
  
  # Methlyation download section --------------------------------------------
  output$download.Methylation.Table = downloadHandler(
    filename = "Methylation.csv",
    content = function(file) {
      write.csv(get.Data.Methylation.Heatmap(input), file,sep=",")
    }
  )
  
  output$download.Plot.Methylation.Heatmap = downloadHandler(
    filename = "Methylation_Heatmap.pdf",
    content = function(file) {
      pdf(file)
      mt = plot.Methylation.Heatmap(input)
      print(mt)
      dev.off()
    }
  )
  
  
  
  
  output$download.Plot.Methylation.Distribution = downloadHandler(
    filename = "Methylation_Distribution_by_Gene.pdf",
    content = function(file) {
      pdf(file)
      plot_output_list = lapply(1:length(input$methylation.Genes), function(i) {
        plotname = paste("plot", i, sep="")
        print(plot.Methylation.Distribution(input,i))
      })
      dev.off()
    }
  )
  
  output$download.Plot.Methylation.Waterfall = downloadHandler(
    filename = "Methylation_Waterfall_by_Gene.pdf",
    content = function(file) {
      pdf(file)
      plot_output_list = lapply(1:length(input$methylation.Genes), function(i) {
        plotname = paste("plot", i, sep="")
        print(plot.Methylation.Waterfall(input,i))
      })
      dev.off()
    }
  )
  
  
  # Expression Section ------------------------------------------------------
  
  # Global Expression -------------------------------------------------------
  
  
  output$plot.Global.Expresion = renderPlot({
    withProgress(session, {
      setProgress(message = "Creating Plot",
                  detail = "This may take a few seconds...")
      mt = plot.Global.Expresion(input)
      tryCatch(print(mt), error=function(e) e)
    })
  })
  
  
  # Boxplots ----------------------------------------------------------------
  
  output$plot.Expression.By.Gene = renderUI({
    number_of_plots = length(input$expression.Genes)
    plot_output_list = lapply(1:number_of_plots, function(i) {
      plotname = paste("plot.Expression.By.Gene", i, sep="")
      plotOutput(plotname, height = hw, width = hw)
    })
    do.call(tagList, plot_output_list)
  })
  
  for (i in 1:max_plots) {
    local({
      my_i = i
      plotname = paste("plot.Expression.By.Gene", my_i, sep="")
      output[[plotname]] = renderPlot({
        withProgress(session, {
          setProgress(message = "Creating Plots",
                      detail = "This may take a few seconds...")
          print(plot.Expression.By.Gene(input, my_i))
        })
      })
    })
  }
  
  output$plot.Expression.by.Entity = renderUI({
    number_of_plots = length(input$expression.Entities)
    plot_output_list = lapply(1:number_of_plots, function(i) {
      plotname = paste("boxplot_by_entity", i, sep="")
      plotOutput(plotname, height = hw, width = hw)
    })
    
    do.call(tagList, plot_output_list)
  })
  
  for (i in 1:max_plots) {
    local({
      my_i = i
      plotname = paste("boxplot_by_entity", my_i, sep="")
      output[[plotname]] = renderPlot({
        withProgress(session, {
          setProgress(message = "Creating Plots",
                      detail = "This may take a few seconds...")
          print(plot.Expression.By.Entity(input, my_i))
        })
      })
    })
  }
  
  
  # Waterfall Plots ---------------------------------------------------------
  
  output$plot.Expression.Waterfall.by.Entity = renderUI({
    number_of_plots = length(input$expression.Entities)
    plot_output_list = lapply(1:number_of_plots, function(i) {
      plotname = paste("plot.Expression.Waterfall.by.Entity", i, sep="")
      plotOutput(plotname, height = hw, width = hw)
    })
    
    do.call(tagList, plot_output_list)
  })
  
  for (i in 1:max_plots) {
    local({
      my_i = i
      plotname = paste("plot.Expression.Waterfall.by.Entity", my_i, sep="")
      output[[plotname]] = renderPlot({
        withProgress(session, {
          setProgress(message = "Creating Plots",
                      detail = "This may take a few seconds...")
          print(plot.Expression.Waterfall.by.Entity(input, my_i))
        })
      })
    })
  }
  
  output$plot.Expression.Waterfall.by.Gene = renderUI({
    number_of_plots = length(input$expression.Genes)
    plot_output_list = lapply(1:number_of_plots, function(i) {
      plotname = paste("plot.Expression.Waterfall.by.Gene", i, sep="")
      plotOutput(plotname, height = hw, width = hw)
    })
    
    do.call(tagList, plot_output_list)
  })
  
  for (i in 1:max_plots) {
    local({
      my_i = i
      plotname = paste("plot.Expression.Waterfall.by.Gene", my_i, sep="")
      output[[plotname]] = renderPlot({
        withProgress(session, {
          setProgress(message = "Creating Plots",
                      detail = "This may take a few seconds...")
          print(plot.Expression.Waterfall.by.Gene(input, my_i))
        })
      })
    })
  }
  
  # Expression CSV ----------------------------------------------------------
  
  output$Expression.Table = renderDataTable({
    genes =  input$expression.Genes
    entities = input$expression.Entities
    ret = generate.Plot.Matrix(genes, entities)
    ret = ret[, c(1:4, ncol(ret))]
    rownames(ret) = NULL
    ret
  }) 
  
  # Expression Download Seciton ---------------------------------------------
  
  output$download.plot.Global.Expresion=downloadHandler(
    filename = "Global_Expression.pdf",
    content = function(file) {
      pdf(file)
      mt = plot.Global.Expresion(input)
      print(mt)
      dev.off()
    }
  )
  
  output$download.plot.Expression.By.Gene=downloadHandler(
    filename = "Expression_by_Gene.pdf",
    content = function(file) {
      pdf(file)
      plot_output_list = lapply(1:length(input$expression.Genes), function(i) {
        plotname = paste("plot", i, sep="")
        print(plot.Expression.By.Gene(input,i))
      })
      dev.off()
    }
  )
  
  output$download.plot.Expression.by.Entity = downloadHandler(
    filename = "Expression_by_Entity.pdf",
    content = function(file) {
      pdf(file)
      plot_output_list = lapply(1:length(input$expression.Entities), function(i) {
        plotname = paste("plot", i, sep="")
        print(plot.Expression.By.Entity(input,i))
      })
      dev.off()
    }
  )
  
  output$download.plot.Expression.Waterfall.by.Entity=downloadHandler(
    filename = "Waterfall_by_Entity.pdf",
    content = function(file) {
      pdf(file)
      plot_output_list = lapply(1:length(input$expression.Entities), function(i) {
        plotname = paste("plot", i, sep="")
        print(plot.Expression.Waterfall.by.Entity(input,i))
      })
      dev.off()
    }
  )
  
  output$download.plot.Expression.Waterfall.by.Gene = downloadHandler(
    filename = "Waterfall_by_Gene.pdf",
    content = function(file) {
      pdf(file)
      plot_output_list = lapply(1:length(input$expression.Genes), function(i) {
        plotname = paste("plot", i, sep="")
        print(plot.Expression.Waterfall.by.Gene(input,i))
      })
      dev.off()
    }
  )
  
  
  output$download.Expression.Table = downloadHandler(
    filename = "Expression.csv",
    content = function(file) {
      genes =  input$expression.Genes
      entities = input$expression.Entities
      ret = generate.Plot.Matrix(genes, entities)
      ret = ret[, c(1:4, ncol(ret))]
      rownames(ret) = NULL
      write.csv(ret, file,sep=",")
    }
  )
  
  
  # Copynumber Section ------------------------------------------------------
  
  output$plot.Global.Copynumber = renderPlot({
    withProgress(session, {
      setProgress(message = "Creating Plot",
                  detail = "This may take a few seconds...")
      mt = plot.Global.Copynumber(input)
      tryCatch(print(mt), error=function(e) e)
    })
  })
  
  
  output$plot.Copynumber.Distribution = renderUI({
    number_of_plots = length(input$copynumber.Entities)
    plot_output_list = lapply(1:number_of_plots, function(i) {
      plotname = paste("plot.Copynumber.Distribution", i, sep="")
      plotOutput(plotname, height = hw, width = hw)
    })
    
    do.call(tagList, plot_output_list)
  })
  
  for (i in 1:max_plots) {
    local({
      my_i = i
      plotname = paste("plot.Copynumber.Distribution", my_i, sep="")
      output[[plotname]] = renderPlot({
        withProgress(session, {
          setProgress(message = "Creating Plots",
                      detail = "This may take a few seconds...")
          print(plot.Copynumber.Distribution(input, my_i))
        })
      })
    })
  }
  
  output$copynumber.Table = renderDataTable({
    genes =  input$copynumber.Genes
    entities = input$copynumber.Entities
    ret = generate.Copynumber.Matrix(genes, entities)
    ret = cbind(ret, Patient=substr(rownames(ret), 1, 12))
    rownames(ret) = NULL
    colnames(ret)[1] = "CNV.Status"
    ret = ret[, c("CNV.Status","Gene","Entity","Patient")]
    ret
  }) 
  
  # CNV Download Section ----------------------------------------------------
  
  output$download.copynumber.Table = downloadHandler(
    filename = "Copynumber.csv",
    content = function(file) {
      genes =  input$copynumber.Genes
      entities = input$copynumber.Entities
      ret = generate.Copynumber.Matrix(genes, entities)
      ret = cbind(ret, Patient=substr(rownames(ret), 1, 12))
      rownames(ret) = NULL
      colnames(ret)[1] = "CNV.Status"
      ret = ret[, c("CNV.Status","Gene","Entity","Patient")]
      write.csv(ret, file,sep=",")
    }
  )
  
  output$plot.Copynumber.Distribution = renderUI({
    number_of_plots = length(input$copynumber.Entities)
    plot_output_list = lapply(1:number_of_plots, function(i) {
      plotname = paste("plot.Copynumber.Distribution", i, sep="")
      plotOutput(plotname, height = hw, width = hw)
    })
    
    do.call(tagList, plot_output_list)
  })
  
  output$download.Plot.Copynumber.Distribution = downloadHandler(
    filename = "Barplot_by_Gene.pdf",
    content = function(file) {
      pdf(file)
      plot_output_list = lapply(1:length(input$copynumber.Entities), function(i) {
        print(plot.Copynumber.Distribution(input,i))
      })
      dev.off()
    }
  )
  
  output$download.Plot.Global.Copynumber = downloadHandler(
    filename = "Heatmap_by_Gene.pdf",
    content = function(file) {
      pdf(file)
      print(plot.Global.Copynumber(input))
      dev.off()
    }
  )  
  
  # Oberserver to Update Views ----------------------------------------------
  
  updateSelectizeInput(session, 'variant.Genes', choices = unique(variant.Table$Hugo_Symbol), server = TRUE)
  updateSelectizeInput(session, 'methylation.Genes', choices = unique(pos.2.gene$UCSC_RefGene_Name), server = TRUE)
  updateSelectizeInput(session, 'expression.Genes', choices = unique(unlist(expression.Table[, "rn" , with=F])), server = TRUE)
  updateSelectizeInput(session, 'copynumber.Genes', choices = unique(unlist(copynumber.Table[, "rn" , with=F])), server = TRUE)
})