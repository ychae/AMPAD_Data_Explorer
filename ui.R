
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
  tags$head(
    singleton(
      includeScript("www/readCookie.js")
    )
  ),
  
  fluidRow(
    column(width = 9,
           # Sample filtering
           fluidRow(height=3,
                    column(width = 12,
                           conditionalPanel('input.custom_search == "Gene_List"',
                           box(width=NULL, solidHeader=TRUE, status="primary",
                               title = tagList(shiny::icon("filter", lib = "glyphicon"), "Filter samples"),
                               collapsible = T,
                               uiOutput("metadataui"),
                               uiOutput("loggedin")
                           )
                    )
                    )      
           ),
           
           # Select your favorite gene input
           box(width = 5, solidHeader = T, status = 'primary', collapsible = T,
               title = tagList(shiny::icon('search', lib = 'glyphicon'), 'Gene Input'),
               selectInput("gene", label = "Select or enter the HGNC symbol of your favorite gene.", 
                           choices = unique(final_data_df$hgnc_symbol), selected = "TSPAN6")
               
           ),
           
           # Redirect to the AMP-AD Knowledge Portal
           box(width = 5, solidHeader = T, status="primary", collapsible = T,
               title = tagList(shiny::icon("info-sign", lib="glyphicon"), 'More Info'),
               a("Click here to go the AMP-AD Knowledge Portal Homepage", 
                 href="https://www.synapse.org/#!Synapse:syn2580853/wiki/66722",
                 target="_blank"))
    ),
           
    
    # Tabs to switch between violin plot, effect size plot and data table
    fluidRow(
      tabBox(
        title = NULL,
        width = NULL, 
        tabPanel('Heatmap',
          box(width = NULL, solidHeader = TRUE,
            plotOutput("heatmap_gene", height = 650)
        )
      ),
        tabPanel("Violin Plots",
                 box( width = NULL, solidHeader = F,
                      fluidRow(
                        splitLayout(cellWidths=c('28.5%', '57.3%', '14.2%'),
                                    plotOutput("violin_plot_mayo"),
                                    plotOutput("violin_plot_msbb"),
                                    plotOutput("violin_plot_rosmap")
                        )),
                      
                      brainRegion, # legend for Tissue abbreviations
                      br(),
                      diseaseStatus,
                      br(),# legend for Status abbreviations
                      downloadButton('downloadPlot')
                 )
        ),
        tabPanel("Effect Size Plots",
                 box( width = NULL, solidHeader = F,
                      fluidRow(
                        
                        plotOutput("forest_plot")
                      )
                 )
        ),                     
        tabPanel("Data Table",
                 dataTableOutput('data_by_gene'))
      )
    ),
#            # Main plot area
#            box(width = NULL, solidHeader = TRUE,
# #                conditionalPanel("input.show_dt",
# #                                 DT::dataTableOutput('infotbl')
# #                ),
# #                conditionalPanel("!input.show_dt",
#                                 plotOutput("heatmap_gene", height = 650)
#            )
#            
#     ),
    
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





