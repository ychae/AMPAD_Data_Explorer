library(DT)
library(ComplexHeatmap)

sample_gene_list <- c("APOE", "CD33", "CD44", "APP", "PER1", 
                      "PICALM", "VGF", "MAPT", "BIN1", "CD47", 
                      "TYROBP", "DOCK2", "FCER1G", "FYN", "MEF2C",
                      "NPTX2", "SNAP29", "STX4", "SYT1", "TARDBP")


#Define the server the logic
shinyServer(
  
  function(input, output, session) {
    
    dataset <- reactive({
      eset.mRNA
    })
    
    output$plotHelp <- renderUI({
      filter_type_text <- filter_type_help()
      p(class = "text-muted", filter_type_text)
    })
    
    output$featureui <- renderUI ({
      featuresel <- input$custom_search
      
      switch(featuresel,
             Gene_List = tagList(p(class = "text-info",
                                  "Enter HGNC gene symbols (e.g., TP53)."),
                                tags$textarea(paste0(c(sample_gene_list), collapse="\n"),
                                              rows=5, id="custom_input_list", style="width: 100%"),
                                actionButton("refreshGene", "Refresh")),
             Significant_Genes = tagList(p(class = "text-info",
                                  "Adjust for desired pvalues and log fold changes"),
                                sliderInput('adjPVal', label=h6('FDR Adjusted p-value'), sep="",
                                            min=0.000001, max=0.05, value=0.05, step=0.001),
                                sliderInput('logFC', label=h6('Log Fold Change'),
                                            min=0, max=4, value=1, step=0.1),
                                actionButton("refreshValue", "Refresh"),
                                hr(),
                                radioButtons("plotdisplay",
                                             label="Select study to plot",
                                             choices=unique(logFC$DataSetName),
                                             selected="Mayo"))
                                
             
    )
    })
    
#     output$plotdisplayui <- renderUI({
#       
#       featuresel <- input$custom_search
#       flog.debug(paste("featuresel = ", featuresel), name='server')
#       
#       if (featuresel %in% c("Gene_List", "Significant_Genes")) {
#         shortfeaturesel <- eset.mRNA
#       }
#       else {
#         shortfeaturesel <- featuresel
#       }
#       
#       radioButtons("plotdisplay",
#                    label="Select study to plot",
#                    choices=c("Mayo", "MSBB", "ROSMAP"),
#                    selected="Mayo")
#     })





filter_type_help <- reactive({
  "Plotting selected genes."
})  

filtered_dataset <- reactive({
  
  adjPVal <- isolate(input$adjPVal)
  logFCThreshold <- isolate(input$logFC)
  ds <- dataset()
  ds_filtered <- filter_by_metadata(input, ds)
  flog.debug(sprintf("filtered ds dims: %s", dim(ds_filtered)), name="server")
  user_feats <- user_submitted_features()
  
  filtered_by_pvalue <- ad_data_pvalue %>% filter(adj.P.Val <= adjPVal, abs(logFC) >= logFCThreshold) 
  
  if (length(input$DataSetName) > 0) {
    filtered_by_pvalue <- filtered_by_pvalue %>% filter(DataSetName %in% input$DataSetName)
  }
  ds_filtered <- ds_filtered[unique(filtered_by_pvalue$ensembl_gene_id),]
  
  
  if (length(user_feats) > 0) {
    userFeats <- intersect(user_feats, fData(ds_filtered)$hgnc_symbol)
    feats <- filter(fData(ds_filtered), hgnc_symbol %in% userFeats)$ensembl_gene_id
    ds_filtered <- ds_filtered[feats, ]
    flog.debug(sprintf("# features in common: %s", length(feats)), name="server")
  }
  
  ds_filtered
})

output$infotbl <- DT::renderDataTable({
  ds <- filtered_dataset()
  rownames(ds) <- fData(ds)$hgnc_symbol
  foo <- signif(exprs(ds), 3)
  # foo <- cbind(feature=featureNames(ds), foo)
  DT::datatable(foo,
                options = list(
                  dom = 'tp',
                  lengthChange = FALSE,
                  pageLength = 15,
                  scrollX = TRUE,
                  scrollCollapse = TRUE))
})


# prepare data for download
output$download_data <- downloadHandler(
  filename = function() {'AMP_AD_data.csv'},
  content  = function(file){
    res <- filtered_dataset()
    hm <- heatmap_cache$heatmap
    
    mat <- exprs(res)[hm$tree_row$order, hm$tree_col$order]        
    
    df <- cbind(data.frame(ID=rownames(mat)),
                as.data.frame(mat))
    write.csv(df, file, row.names=F, col.names=T)
  }
)

user_submitted_features <- reactive({
  if (input$custom_search == "Gene_List") {
    input$refreshGene
  }
  else if(input$custom_search == "Significant_Genes") {
    input$refreshValue
  }
  
  geneList <- isolate(input$custom_input_list)
  
  curr_filter_type <- paste(input$custom_search, input$plotdisplay, sep="_")
  flog.debug(curr_filter_type, name="server")
  
  featureList <- clean_list(geneList, change_case=toupper)
  
  flog.debug(sprintf("In %s, selected %s features", curr_filter_type, length(featureList)), name="server")
  
  featureList
})

#     output$featxsamples <- renderInfoBox({
#       ds <- filtered_dataset()
#       infoBox(title="Features x Samples", 
#               value=sprintf("%s x %s", nrow(ds), ncol(ds)),
#               fill=TRUE, width=NULL)
#     })

heatmap_cache <- reactiveValues()

#return the heatmap plot
output$heatmap <- renderPlot({  
  flog.debug("Making heatmap", name='server')
  
  m_eset <- filtered_dataset()
  m <- exprs(m_eset)
  m <- data.matrix(m)
  m <- m[, -c(2:6, 8:12, 17, 19)] #Keep only the AD vs Controls
  
  rownames(m) <- fData(m_eset)$hgnc_symbol
  
  validate( need( ncol(m) != 0, "Filtered matrix contains 0 Samples.") )
  validate( need( nrow(m) != 0, "Filtered matrix contains 0 features.") )
  validate( need(nrow(m) < 15000, "Filtered matrix contains > 10000 genes.") )
  
  filtered_metadata <- pData(m_eset)
  annotation <- get_heatmapAnnotation(input$heatmap_annotation_labels, filtered_metadata)
  
  fontsize_row <- ifelse(nrow(m) > 100, 0, 8)
  fontsize_col <- ifelse(ncol(m) > 50, 0, 8)    
  
  # Need to scale methylation breaks differently
  heatmap.color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  ha <- HeatmapAnnotation(df = annotation)
  print(ha)
  h <- memoise(Heatmap(matrix = m, 
                       cluster_rows = FALSE, 
                       cluster_columns = FALSE, 
                       col = heatmap.color,
                       heatmap_legend_param = list(title = "Log Fold \nChange", 
                                                   title_position = "topcenter",
                                                   ncol = 1),
                       top_annotation = ha
  ))
  heatmap_cache$heatmap <- print(h)
  
})
  }
)


