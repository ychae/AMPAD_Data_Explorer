
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


myHeader <- dashboardHeader(title="AMP-AD Data Explorer", disable=T)

mySidebar <- dashboardSidebar(disable=TRUE)

myBody <-dashboardBody(
  fluidRow(
    column(width = 9,
           
           # Sample filtering
           fluidRow(height=3,
                    column(width = 12,
                           conditionalPanel('input.custom_search == "Gene_List"',
                           box(width=NULL, solidHeader=TRUE, status="primary",
                               title = tagList(shiny::icon("filter", lib = "glyphicon"), "Filter samples"),
                               collapsible = T,
                               tags$table(class="table table-condensed",
                                          tags$tr(
                                            tags$td(selectInput('BrainRegion', h6('Brain Region'),
                                                                choices=unique(covariates$BrainRegion),
                                                                selectize=T, multiple=T)),
                                            tags$td(selectInput('Study', h6('Study'),
                                                                choices=unique(covariates$Study),
                                                                selectize=T, multiple=T))
                                          ),
                                          tags$tr(

                                            tags$td(selectInput('Status', h6('Diagnosis'),
                                                                choices=unique(covariates$Status),
                                                                selectize=T, multiple=T)),
                                            tags$td(selectInput('Gender', h6('Gender'),
                                                                choices=unique(covariates$Gender),
                                                                selectize=T, multiple=T))
                                            )
                               
                               )
                           )
                    )
                    )      
           ),
           
           # Main plot area
           box(width = NULL, solidHeader = TRUE,
#                conditionalPanel("input.show_dt",
#                                 DT::dataTableOutput('infotbl')
#                ),
#                conditionalPanel("!input.show_dt",
                                plotOutput("heatmap_gene", height = 650)
           )
           
    ),
    
    column(width = 3,
           
           # Plot selection box
           box(width = NULL, status = "primary", solidHeader=TRUE,
               title = "Select features to display", collapsible = T,
               
               selectInput("custom_search",
                           label = "Select feature type",
                           choices= c("Gene_List", "Significant_Genes"),
                           selectize = T, multiple = F, selected = "Gene_List"),
               
               uiOutput("featureui")
               
               #hr(),
               
               #checkboxInput('show_dt', 'Show data values instead of heatmap', value = FALSE)
               #uiOutput("plotHelp")               
           )
           
           
           # Download box
#            box(width=NULL, status = 'info', solidHeader=TRUE,
#                collapsible=TRUE, collapsed=FALSE,
#                title = tagList(shiny::icon("save", lib = "glyphicon"), "Download data"),
#                selectInput("savetype",
#                            label=h6("Save as:"),
#                            choices=c("comma separated (CSV)", "tab separated (TSV)"),
#                            selectize=F, multiple=F, selected="comma separated (CSV)"),
#                downloadButton(outputId='download_data', label='Download')
#            )
    )
  )
)

dashboardPage(header=myHeader, sidebar=mySidebar, body=myBody,
              skin = "blue")





