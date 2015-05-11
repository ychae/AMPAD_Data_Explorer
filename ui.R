# library(shinyIncubator)

meth_data_notes <- '<pre>Data Processing Notes:<br>Methylation probes with variation &gt; .01 across all samples were choosen from the normalized data matrix(<a href="https://www.synapse.org/#!Synapse:syn2233188" target="_blank">syn223318</a>). The probes were selected based on genes using a mapping file.(<a href="https://www.synapse.org/#!Synapse:syn2324928" target="_blank">syn2324928</span></a>). Hierarchical clustering was used to cluster rows and columns.</pre>'

#2. mRNA data notes
mRNA_data_notes  <- 'Data Processing Notes:<br>Using mRNA normalized data matrix from <a href="https://www.synapse.org/#!Synapse:syn2701943" target="_blank">syn2701943</a> and metadata from <a href="https://www.synapse.org/#!Synapse:syn2731147" target="_blank">syn2731147</a>. Hierarchical clustering was used to cluster rows and columns.'

#3. miRNA data notes
miRNA_data_notes <- 'Data Processing Notes:<br>Using miRNA normalized data matrix from <a href="https://www.synapse.org/#!Synapse:syn2701942" target="_blank">syn2701942</a> and metadata from <a href="https://www.synapse.org/#!Synapse:syn2731149" target="_blank">syn2731149</a>. The miRNAs were selected based on target genes using a mapping file <a href="https://www.synapse.org/#!Synapse:syn2246991" target="_blank">syn2246991</a>. Hierarchical clustering was used to cluster rows and columns.'

#main UI code

# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shinydashboard)

myHeader <- dashboardHeader(title="PCBC Data Explorer", disable=TRUE)

mySidebar <- dashboardSidebar(disable=TRUE)
#   sidebarMenu(
#     menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard"))
#   ))

myBody <-dashboardBody(
  fluidRow(
    column(width = 9,
           fluidRow(height=3,
                    column(width = 9,
                           box(width=NULL, solidHeader=TRUE, status="primary",
                               title="Sample Filters",
                               tags$table(class="table table-condensed",
                                          tags$tr(
                                            tags$td(selectInput('linetype', h6('Line type'),
                                                                choices=c("Bzzzz", "Fizzzzz"), #unique(combined_metadata$Cell_Line_Type),
                                                                selectize=T, multiple=T, selected=c('ESC','iPSC'))),
                                            tags$td(selectInput('vector_type', h6('Reprogramming Vector'),
                                                                choices=c("Baz", "Zoo"), #unique(combined_metadata$Reprogramming_Vector_Type),
                                                                selectize=T, multiple=T)),
                                            tags$td(selectInput('gene_combination', h6('Reprogramming Genes'),
                                                                choices=c("Boo", "Far"), #unique(combined_metadata$Reprogramming_Gene_Combination),
                                                                selectize=T, multiple=T)),
                                            tags$td(selectInput('tissue_origin', h6('Tissue'),
                                                                choices=c("Foo", "Bar"),#unique(combined_metadata$Tissue_of_Origin),
                                                                selectize=T, multiple=T))
                                          ),
                                          
                                          tags$tr(
                                            tags$td(selectInput('diff_state', h6('Differentiation'),
                                                                choices=c("Foo", "Bar"), #unique(combined_metadata$Diffname_short),
                                                                selectize=T, multiple=T)),
                                            tags$td(selectInput('cell_origin', h6('Cell Type of Origin'),
                                                                choices=c("Foo", "Bar"), #unique(combined_metadata$Cell_Type_of_Origin),
                                                                selectize=T, multiple=T)),
                                            tags$td(selectInput('gender', h6('Gender'),
                                                                choices=c("Foo", "Bar"), #unique(combined_metadata$Gender),
                                                                selectize=T, multiple=T)),
                                            tags$td(selectInput('originating_lab', h6('Originating Lab'),
                                                                choices=c("Foo", "Bar", "Some really long lab name with many people"), #unique(combined_metadata$Originating_Lab_ID),
                                                                selectize=T, multiple=T))
                                          )           
                               )
                           )
                    ),
                    column(width = 3,
                           
                           box(width=NULL, status='primary', collapsible=TRUE, solidHeader=TRUE,
                               title="Sample labeling",
                               selectInput('heatmap_annotation_labels',
                                           'Annotate Samples by:',
                                           choices=c("Diffname_short", "Foo", "Bar"), #colnames(combined_metadata)[-1],  #-1 to remove the first value "Sample"
                                           selected='Diffname_short')               
                           )
                    )                       
           ),
           box(width = NULL, solidHeader = TRUE,
               plotOutput("plot1", height = 700))
           
    ),
    #)
    
    column(width = 3,
           box(width = NULL, status = "primary", solidHeader=TRUE,
               title="Plot to display",
               selectInput("plotdisplay",
                           label=NULL, #h6(""),
                           choices=c("mRNA", "miRNA", "Methylation"),
                           selectize=T, multiple=F, selected="mRNA"),
               
               conditionalPanel('(input.search_box == "Gene" | input.search_box == "Pathway") & input.plotdisplay == "mRNA"',
                                p("Plotting selected genes.")),
               
               conditionalPanel('(input.search_box == "Gene" | input.search_box == "Pathway") & input.plotdisplay == "miRNA"',
                                p("Plotting miRNAs targeting selected genes.")),
               
               conditionalPanel('(input.search_box == "Gene" | input.search_box == "Pathway") & input.plotdisplay == "Methylation"',
                                p("Plotting methylation probes of selected genes.")),
               
               conditionalPanel('input.search_box == "miRNA" & input.plotdisplay == "miRNA"',
                                p("Plotting selected miRNAs.")),
               
               conditionalPanel('input.search_box == "miRNA" & input.plotdisplay == "mRNA"',
                                p("Plotting genes targeted by selected miRNAs.")),
               
               conditionalPanel('input.search_box == "miRNA" & input.plotdisplay == "Methylation"',
                                p("Plotting methylation probes for genes targeted by selected miRNAs."))
               
           ),
           
           tabBox(width=NULL, status="info",
                  id="search_box",
                  # Title can include an icon
                  title = tagList(shiny::icon("search")),
                  tabPanel("Gene",
                           tags$textarea(paste0(sample_gene_list, collapse="\n"),
                                         rows=5, id="Genes", style="width: 100%"),
                           p(class = "text-muted",
                             "Accepts HUGO/Ensembl/Entrez gene ids")),
                  tabPanel("Pathway", 
                           selectInput("selected_pathways", label=NULL,
                                       choices = c("Pathway 1", "Pathway 2"), #names(pathways_list),
                                       selectize=F, multiple=F)),
                  tabPanel("miRNA", 
                           tags$textarea(paste0(sample_mirna_list, collapse="\n"),
                                         rows=5, id="miRNA", style="width: 100%")),
                  actionButton("Refresh", "Refresh"),
           ),
           box(width = NULL, status = "warning", solidHeader=TRUE, collapsible=TRUE,
               title="Correlation",
               
               conditionalPanel('input.search_box == "miRNA"',
                                "Not available."),
               
               conditionalPanel('input.search_box != "miRNA"',
                                checkboxInput('incl_corr_genes', 
                                              'also include correlated genes', 
                                              value = FALSE),
                                
                                conditionalPanel(
                                  condition="input.incl_corr_genes",
                                  sliderInput('corr_threshold', label=h6('Correlation Threshold'),
                                              min=0.5, max=1.0, value=0.9, step=0.05),
                                  # correlation direction
                                  selectInput("correlation_direction",
                                              label=h6("Correlation Direction"),
                                              choices=c("both", "positive", "negative"),
                                              selectize=T, multiple=F, selected="both"),
                                  p(class = "text-muted",
                                    br(),
                                    "This is an example note in a muted text color."
                                  )
                                )
               )
           ),
           box(width = NULL, status = "warning", collapsible=TRUE, solidHeader=TRUE,
               title="Clustering",
               #distance metric
               selectInput("clustering_distance", "Distance Calculation",
                           choices=c("correlation", "euclidean", "maximum", 
                                     "manhattan", "canberra", "binary", "minkowski"),
                           selectize=T, multiple=F, selected="correlation"),
               
               #linkage 
               selectInput("clustering_method", "Clustering Method",
                           choices=c("ward", "single", "complete", "average", 
                                     "mcquitty", "median", "centroid"),
                           selectize=T, multiple=F, selected="complete"),
               
               checkboxInput('cluster_cols', 'Cluster the columns', value = TRUE),
               
               checkboxInput('cluster_rows', 'Cluster the rows', value = TRUE)
               
           ),
           
           box(width=NULL, status = 'info', collapsible=TRUE, solidHeader=TRUE,
               title = tagList(shiny::icon("save", lib = "glyphicon"), "Download"),
               selectInput("savetype",
                           label=h6("Save as:"),
                           choices=c("comma separated (CSV)", "tab separated (TSV)"),
                           selectize=F, multiple=F, selected="comma separated (CSV)"),
               actionButton("Download", "Download")
           )
    )
  )
)

# myBody <- dashboardBody(
#   tabItems(
#     # First tab content
#     tabItem(tabName = "dashboard",
#             fluidRow(
#               box(
#                 title = "Filter Features", width = 2, solidHeader = TRUE,
#                 )),
#             fluidRow(
#               box(
#                 title = "Filter Samples", width = 10, solidHeader = TRUE, status = "primary",
#                 "Box content"
#               ))
#             )
#   )
# )

dashboardPage(header=myHeader, sidebar=mySidebar, body=myBody,
              skin = "blue")


# tags$textarea(id="custom_gene_list",
#               rows=8, cols=50,
#               paste0(sample_gene_list, collapse=', ')),
# 
# 
# h5('1.b. Add miRNA Targets (mirbase ids):'),
# tags$textarea(id="custom_miRNA_list",rows=4,cols=50),                 
# 
# actionButton("custom_search", h4("Update")),
# value='custom_gene_list'
#         
# # TAB PANEL 2 : select a pathway
# selectInput("selected_pathways",
#             h5("1.a. Select Pathway/s"),
#             choices = names(pathways_list),
#             selectize=T, multiple=T, width='400px',
#             selected = names(pathways_list)[c(1:2)])
#       
# #Main shiny panel
# plotOutput("mRNA_heatMap",height="700px",width="auto",hoverId=NULL),
# htmlOutput("topgene_linkOut"),
# downloadButton('download_mRNAData','Download mRNA expression data'),
# HTML(mRNA_data_notes)
# plotOutput("microRNA_heatMap",height="700px",width="auto",hoverId=NULL),
# downloadButton('download_miRNAData','Download microRNA expression data'),
# HTML(miRNA_data_notes)
# 
# plotOutput("methylation_heatMap",height="700px",width="auto",hoverId=NULL),
# downloadButton('download_methylationData','Download methylation data'),
# HTML(meth_data_notes)
