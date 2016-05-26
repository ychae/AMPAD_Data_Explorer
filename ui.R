
#meth_data_notes <- '<pre>Data Processing Notes:<br>Methylation probes with variation &gt; .01 across all samples were choosen from the normalized data matrix(<a href="https://www.synapse.org/#!Synapse:syn2233188" target="_blank">syn223318</a>). The probes were selected based on genes using a mapping file.(<a href="https://www.synapse.org/#!Synapse:syn2324928" target="_blank">syn2324928</span></a>). Hierarchical clustering was used to cluster rows and columns.</pre>'

#2. mRNA data notes
mRNA_data_notes  <- 'Data Processing Notes:<br>Using mRNA normalized data matrix from <a href="https://www.synapse.org/#!Synapse:syn2701943" target="_blank">syn2701943</a> and metadata from <a href="https://www.synapse.org/#!Synapse:syn2731147" target="_blank">syn2731147</a>. Hierarchical clustering was used to cluster rows and columns.'

#3. miRNA data notes
#miRNA_data_notes <- 'Data Processing Notes:<br>Using miRNA normalized data matrix from <a href="https://www.synapse.org/#!Synapse:syn2701942" target="_blank">syn2701942</a> and metadata from <a href="https://www.synapse.org/#!Synapse:syn2731149" target="_blank">syn2731149</a>. The miRNAs were selected based on target genes using a mapping file <a href="https://www.synapse.org/#!Synapse:syn2246991" target="_blank">syn2246991</a>. Hierarchical clustering was used to cluster rows and columns.'

#main UI code
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shinydashboard)
library(DT)

# sample_gene_list <- c("APOE", "CD33", "CD44", "APP", "PER1", 
#                       "PICALM", "VGF", "MAPT", "BIN1", "CD47", 
#                       "TYROBP", "DOCK2", "FCER1G", "FYN", "MEF2C",
#                       "NPTX2", "SNAP29", "STX4", "SYT1", "TARDBP")

myHeader <- dashboardHeader(title="AMP-AD Data Explorer", disable=T)

mySidebar <- dashboardSidebar(disable=TRUE)

myBody <-dashboardBody(
  fluidRow(
    column(width = 9,
           
           # Sample filtering
           fluidRow(height=3,
                    column(width = 12,
                           box(width=NULL, solidHeader=TRUE, status="primary",
                               title = tagList(shiny::icon("filter", lib = "glyphicon"), "Filter samples"),
                               tags$table(class="table table-condensed",
                                          tags$tr(
                                            tags$td(selectInput('DataSetName', h6('Study'),
                                                                choices=unique(logFC$DataSetName),
                                                                selectize=T, multiple=T))
                                          )           
                               )
                           )
                    )
                    
           ),
           
           # Main plot area
           box(width = NULL, solidHeader = TRUE,
               conditionalPanel("input.show_dt",
                                DT::dataTableOutput('infotbl')),
               
               conditionalPanel("!input.show_dt",
                                plotOutput("heatmap", height = 650))
           )
           
    ),
    
    column(width = 3,
           
           # Choose sample labels
           box(width=NULL, status='primary', collapsible=TRUE, 
               collapsed=TRUE, solidHeader=TRUE,
               title = tagList(shiny::icon("th-list", lib="glyphicon"),
                               "Label samples"),               
               selectInput('heatmap_annotation_labels',
                           'Annotate Samples by:',
                           # -1 to remove the first value "Sample"
                           choices=colnames(combined_metadata)[-c(1,4)],
                           selected='DataSetName')
           ),
           
           
           # Plot selection box
           box(width = NULL, status = "primary", solidHeader=TRUE,
               title = "Select features to display",
               
               selectInput("custom_search",
                           label = "Select feature type",
                           choices= c("Gene_List", "Significant_Genes"),
                           selectize = T, multiple = F, selected = "Significant_Genes"),
               
               uiOutput("featureui"),
               
               hr(),
               
               
#                tagList(p(class = "text-info",
#                          "Enter HGNC gene symbols (e.g., TP53)."),
#                        tags$textarea(paste0(c(sample_gene_list), collapse="\n"),
#                                      rows=5, id="custom_input_list", style="width: 100%"),
#                        sliderInput('adjPVal', label=h6('FDR Adjusted p-value'), sep="",
#                                    min=0.000001, max=0.05, value=0.05, step=0.001),
#                        sliderInput('logFC', label=h6('Log Fold Change'),
#                                    min=0, max=4, value=1, step=0.1),
#                        actionButton("refreshGene", "Refresh")),
#                
               
              # uiOutput("plotdisplayui"),
               
               #hr(),
               
               checkboxInput('show_dt', 'Show data values instead of heatmap', value = FALSE)
               #uiOutput("plotHelp")               
           ),
           
           
           # Download box
           box(width=NULL, status = 'info', solidHeader=TRUE,
               collapsible=TRUE, collapsed=FALSE,
               title = tagList(shiny::icon("save", lib = "glyphicon"), "Download data"),
               selectInput("savetype",
                           label=h6("Save as:"),
                           choices=c("comma separated (CSV)", "tab separated (TSV)"),
                           selectize=F, multiple=F, selected="comma separated (CSV)"),
               downloadButton(outputId='download_data', label='Download')
           )
    )
  )
)

dashboardPage(header=myHeader, sidebar=mySidebar, body=myBody,
              skin = "blue")



