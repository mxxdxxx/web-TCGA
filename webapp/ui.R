shinyUI(
  navbarPage("Shiny TCGA - choose analyte type",footer=list(div(HTML(paste("Build: ",system("git rev-list HEAD --count", intern=T),sep="")))) ,
             tabPanel("Variant", pageWithSidebar(
               
               headerPanel(""),
               
               sidebarPanel(
                 
                 tabsetPanel(
                   tabPanel("Cancer Types", checkboxGroupInput("variant.Entities","", choices=sort(unique(variant.Table$Entity)))),
                   tabPanel("Genes", selectizeInput("variant.Genes", NULL, NULL, multiple = T, options = NULL)),
                   tabPanel("Options",
                            checkboxInput(inputId="variant.Exclude.Silent",label="Exclude Silent Mutations",value=T),
                            checkboxInput(inputId="variant.Multi.Mut.Per.Gene",label="Use multiple Mutations per Gene",value=T),
                            wellPanel(
                              p(strong("Only used for Pie-Charts and Table")),
                              numericInput(inputId="threshold",label="Treshold for excluding Mutation (%)",value="0"),
                              checkboxInput(inputId="mut_type_split",label="Split Variant by Classifcation",value=F))
                   )
                 ),
                 wellPanel(submitButton("Run Analysis"))
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("Heatmap", plotOutput("plot.Variant.Heatmap"),downloadButton("download.Plot.Variant.Heatmap", "Download Heatmap")),
                   tabPanel("Pie-Charts", uiOutput("plot.Variant.PieCharts"), downloadButton("download.Plot.Variant.PieCharts", "Download Pie-Charts")),
                   tabPanel("CSV Output", dataTableOutput("variant.Table"), downloadButton("download.Variant.Table", "Download Variants"))
                 )
               )
             )),
             tabPanel("Methylation", pageWithSidebar(
               
               headerPanel(""),
               
               sidebarPanel(
                 
                 tabsetPanel(
                   tabPanel("Cancer Entity", selectInput("methylation.Entity","",choices=sort(names(methylation.List)),selectize=F)),
                   tabPanel("Genes", selectizeInput("methylation.Genes",NULL, NULL, multiple = T,options = NULL)),
                   tabPanel("Options",
                            checkboxInput(inputId="methylation.Paired.Only",label="Only consider paired samples for analysis",value=T),
                            checkboxInput(inputId="methylation.Use.Percentage",label="Use percentages instead of hard numbers",value=T),
                            numericInput(inputId="methylation.Threshold",label="Threshold for differential methylation", value="0.2", min=0),
                            selectInput("methylation.P.or.Diff","",choices=c("P-Value", "Difference"),selectize=F)
                   )
                 ),
                 wellPanel(submitButton("Run Analysis"))
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("Global Differential Methylation", plotOutput("plot.Methylation.Heatmap"),downloadButton("download.Plot.Methylation.Heatmap", "Download Map")),
                   tabPanel("Methylation Distribution", uiOutput("plot.Methylation.Distribution"), downloadButton("download.Plot.Methylation.Distribution", "Download Histograms")),
                   tabPanel("Differential Methylated Samples", uiOutput("plot.Methylation.Waterfall"), downloadButton("download.Plot.Methylation.Waterfall", "Download Waterfallplots")),
                   tabPanel("CSV Output", dataTableOutput("methylation.Table"), downloadButton("download.Methylation.Table", "Download CSV"))
                   
                 )
               )
             )),
             tabPanel("Expression", pageWithSidebar(
               
               headerPanel(""),
               
               sidebarPanel(
                 
                 tabsetPanel(
                   tabPanel("Cancer Types", checkboxGroupInput("expression.Entities","",choices=sort(unique(getCancerTypes(colnames(expression.Table)))))),
                   tabPanel("Genes", selectizeInput("expression.Genes",NULL, NULL, multiple = T,options = NULL)),
                   tabPanel("Options",
                            numericInput(inputId="expression.Threshold",label="Threshold +/-", value="2", min=0),
                            selectInput("expression.Z.or.Foldchange","",choices=c("Z-Score", "Foldchange"),selectize=F),
                            checkboxInput(inputId="expression.Use.Percentage",label="Use percentages instead of hard numbers",value=T)
                            )
                 ),
                 wellPanel(submitButton("Run Analysis"))
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("Global Expression", progressInit(), plotOutput("plot.Global.Expresion"), downloadButton("download.plot.Global.Expresion", "Download Map")),
                   tabPanel("Expression by Gene", progressInit(), uiOutput("plot.Expression.By.Gene"), downloadButton("download.plot.Expression.By.Gene", "Download Boxplots")),
                   tabPanel("Expression by Entity", progressInit(), uiOutput("plot.Expression.by.Entity"), downloadButton("download.plot.Expression.by.Entity", "Download Boxplots")),
                   tabPanel("Waterfall by Entity", progressInit(), uiOutput("plot.Expression.Waterfall.by.Entity"), downloadButton("download.plot.Expression.Waterfall.by.Entity", "Download Waterfall")),
                   tabPanel("Waterfall by Gene", progressInit(), uiOutput("plot.Expression.Waterfall.by.Gene"), downloadButton("download.plot.Expression.Waterfall.by.Gene", "Download Waterfall")),
                   tabPanel("CSV Output", dataTableOutput("expression.Table"), downloadButton("download.Expression.Table", "Download CSV"))
                 )
               )
             )),
             tabPanel("Copy Number Variation", pageWithSidebar(
               
               headerPanel(""),
               
               sidebarPanel(
                 
                 tabsetPanel(
                   tabPanel("Cancer Types", checkboxGroupInput("copynumber.Entities","",choices=sort(unique(getCancerTypes(colnames(expression.Table)))))),
                   tabPanel("Genes", selectizeInput("copynumber.Genes",NULL, NULL, multiple = T,options = NULL)),
                   tabPanel("Options",
                            selectInput("copynumber.Hi.Low","Amplifications to Use",choices=c("Both", "Low Level Only", "High Level Only"),selectize=F),
                            checkboxInput(inputId="copynumber.Use.Percentage",label="Use percentages instead of hard numbers",value=T)
                   )
                 ),
                 wellPanel(submitButton("Run Analysis"))
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("Global Copynumber Variation", progressInit(), plotOutput("plot.Global.Copynumber"), downloadButton("download.Plot.Global.Copynumber", "Download Map")),
                   tabPanel("Copynumber Entity Rate", progressInit(), uiOutput("plot.Copynumber.Distribution"), downloadButton("download.Plot.Copynumber.Distribution", "Download Barplots")),
                   tabPanel("CSV Output", dataTableOutput("copynumber.Table"), downloadButton("download.copynumber.Table", "Download CSV"))
                 )
               )
             )),
             tabPanel("Info", fluidPage(
               titlePanel(""),
                                      includeMarkdown("r-objects/Stats.md")
                                      )
                      
              )
  )
)