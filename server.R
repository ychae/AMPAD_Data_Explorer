library(DT)

clean_list <- function(x, change_case=toupper) {
  # Split by space, comma or new lines
  x <- unlist(strsplit(x, split=c('[\\s+,\\n+\\r+)]'),perl=T))
  
  # convert everything to specified case
  x <- change_case(x)
  
  # remove the blank entries
  x <- x[!(x == "")]
  
  x
}

#Define the server the logic
shinyServer(
  
  function(input, output, session) {
    
    dataset <- reactive({
      switch(input$plotdisplay,
             mRNA = eset.mRNA,
             miRNA = eset.miRNA,
             Methylation = eset.meth)
      
    })
    
    output$plotHelp <- renderUI({
      filter_type_text <- filter_type_help()
      p(class = "text-muted", filter_type_text)
    })
    
    filter_type_help <- reactive({
      curr_filter_type <- paste(input$custom_search, input$plotdisplay, sep="_")
      
      switch(curr_filter_type,
             Gene_mRNA="Plotting selected genes.",
             Pathway_mRNA="Plotting selected genes.",
             miRNA_mRNA="Plotting genes targeted by selected miRNAs.",
             
             Gene_miRNA="Plotting miRNAs targeting selected genes.",
             Pathway_miRNA="Plotting miRNAs targeting selected genes.",
             miRNA_miRNA="Plotting selected miRNAs.",
             
             Gene_Methylation="Plotting miRNAs targeting selected genes.",
             Pathway_Methylation="Plotting miRNAs targeting selected genes.",
             miRNA_Methylation="Plotting methylation probes for genes targeted by selected miRNAs.",
             "Unknown selection.")
      
    })  
    
    filtered_dataset <- reactive({
      
      ds <- dataset()
      ds_filtered <- filter_by_metadata(input, ds)
      
      user_feats <- user_submitted_features()
      feats <- intersect(user_feats, featureNames(ds_filtered))
      flog.debug(sprintf("# features in common: %s", length(feats)), name="server")
      ds_filtered <- ds_filtered[feats, ]
      
      if (input$incl_corr_genes == 'TRUE' & input$plotdisplay == 'mRNA' & 
            input$custom_search %in% c("Gene", "Pathway")) { 
        
        ds_filtered <- get_eset_withcorrelated_genes(feats, dataset(),
                                                     input$corr_threshold,
                                                     input$correlation_direction)
      }
      
      # zero variance filter
      rows_to_keep <- apply(exprs(ds_filtered), 1, var) > 0
      ds_filtered <- ds_filtered[rows_to_keep, ]
      
      ds_filtered
    })
    
    output$infotbl <- DT::renderDataTable({
      ds <- filtered_dataset()
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
      filename = function() {'PCBC_data.csv'},
      content  = function(file){
        res <- filtered_dataset()
        mat <- exprs(res)        
        df <- cbind(data.frame(ID=rownames(mat)),
                    as.data.frame(mat))
        write.csv(df, file, row.names=F, col.names=T)
      }
      )
    
    user_submitted_features <- reactive({
      if (input$custom_search == "Gene") {
        input$refreshGene
      }
      else if(input$custom_search == "miRNA") {
        input$refreshmiRNA
      }
      
      geneList <- isolate(input$custom_input_list)
      selectedPathway <- isolate(input$selected_pathways)
      mirnaList <- isolate(input$custom_mirna_list)
      
      curr_filter_type <- paste(input$custom_search, input$plotdisplay, sep="_")
      flog.debug(curr_filter_type, name="server")
      
      if (curr_filter_type == "Gene_mRNA") {
        featureList <- clean_list(geneList, change_case=toupper)
        featureList <- convert_to_ensemblIds(featureList)
      }
      else if (curr_filter_type == "Pathway_mRNA") {
        featureList <- as.character(unlist(pathways_list[selectedPathway]))
        featureList <- clean_list(featureList, change_case=toupper)
        featureList <- convert_to_ensemblIds(featureList)
      }
      else if (curr_filter_type == "Gene_miRNA") {
        featureList <- clean_list(geneList, change_case=toupper)
        featureList <- convert_to_ensemblIds(featureList)
        selected_miRNAs <- filter(miRNA_to_genes, GeneID %in% featureList)
        featureList <- unique(selected_miRNAs$original)
      }
      else if (curr_filter_type == "Pathway_miRNA") {
        featureList <- as.character(unlist(pathways_list[selectedPathway]))
        featureList <- clean_list(featureList, change_case=toupper)
        featureList <- convert_to_ensemblIds(featureList)
        selected_miRNAs <- filter(miRNA_to_genes, GeneID %in% featureList)
        featureList <- unique(selected_miRNAs$original)
      }
      else if(curr_filter_type == "miRNA_miRNA") {
        featureList <- clean_list(mirnaList, change_case=tolower)
        selected_miRNAs <- filter(miRNA_to_genes, miRNAPrecursor %in% featureList | miRNA1 %in% featureList | 
                                    miRNA2 %in% featureList)
        featureList <- unique(selected_miRNAs$original)
      }
      else if(curr_filter_type == "miRNA_mRNA") {
        featureList <- clean_list(mirnaList, change_case=tolower)
        selected_miRNAs <- filter(miRNA_to_genes, miRNAPrecursor %in% featureList | miRNA1 %in% featureList | 
                                    miRNA2 %in% featureList)
        featureList <- unique(selected_miRNAs$GeneID)
      }
      else if (curr_filter_type == "Gene_Methylation") {
        featureList <- clean_list(geneList, change_case=toupper)
        featureList <- convert_to_EntrezIds(featureList)
        flt_res <- filter(meth_to_gene, entrezID %in% featureList)
        featureList <- unique(flt_res$methProbe)
      }
      else if (curr_filter_type == "miRNA_Methylation") {
        featureList <- clean_list(mirnaList, change_case=tolower)
        selected_miRNAs <- filter(miRNA_to_genes, miRNAPrecursor %in% featureList | miRNA1 %in% featureList | 
                                    miRNA2 %in% featureList)
        featureList <- unique(selected_miRNAs$GeneID)
        featureList <- convert_to_EntrezIds(featureList)
        flt_res <- filter(meth_to_gene, entrezID %in% featureList)
        featureList <- unique(flt_res$methProbe)
      }
      else {
        featureList <- c()
      }
      
      flog.debug(sprintf("In %s, selected %s features", curr_filter_type, length(featureList)), name="server")
      
      featureList
    })
   
    output$featxsamples <- renderInfoBox({
      ds <- filtered_dataset()
      infoBox(title="Features x Samples", 
              value=sprintf("%s x %s", nrow(ds), ncol(ds)),
              fill=TRUE, width=NULL)
    })
    
    heatmap_cache <- reactiveValues()
    
    #return the heatmap plot
    output$heatmap <- renderPlot({  
      flog.debug("Making heatmap", name='server')
      
      cluster_rows <- input$cluster_rows
      cluster_cols <- input$cluster_cols
      
      m_eset <- filtered_dataset()
      m <- exprs(m_eset)
      m <- data.matrix(m)
      
      validate( need( ncol(m) != 0, "Filtered mRNA expression matrix contains 0 Samples") )
      validate( need( nrow(m) != 0, "Filtered mRNA expression matrix contains 0 genes") )
      validate( need(nrow(m) < 10000, "Filtered mRNA expression matrix contains > 10000 genes. MAX LIMIT 10,000 ") )
      
      filtered_metadata <- pData(m_eset)
      annotation <- get_heatmapAnnotation(input$heatmap_annotation_labels, filtered_metadata)
      
      fontsize_row <- ifelse(nrow(m) > 100, 0, 8)
      fontsize_col <- ifelse(ncol(m) > 50, 0, 8)    
      
      withProgress(session, {
        setProgress(message = "clustering & rendering heatmap, please wait", 
                    detail = "This may take a few moments...")
        heatmap_cache$heatmap <- expHeatMap(m,annotation,
                                            clustering_distance_rows = input$clustering_distance,
                                            clustering_distance_cols = input$clustering_distance,
                                            fontsize_col=fontsize_col, 
                                            fontsize_row=fontsize_row,
                                            scale=T,
                                            clustering_method = input$clustering_method,
                                            explicit_rownames = fData(m_eset)$explicit_rownames,
                                            cluster_rows=cluster_rows, cluster_cols=cluster_cols,
                                            drawColD=FALSE)
      }) #END withProgress
    })
    
  }
)


#   output$topgene_linkOut <- reactive({
#     prefix <- '<form action="https://toppgene.cchmc.org/CheckInput.action" method="post" target="_blank" display="inline">\
#     <input type="hidden" name="query" value="TOPPFUN">\
#     <input type="hidden" id="type" name="type" value="HGNC">\
#     <input type="hidden" name="training_set" id="training_set" value="%s">\
#     <input type="Submit" class="btn shiny-download-link" value="Enrichment Analysis in ToppGene">\
#     </form>'
#     geneIds <- rownames(get_filtered_mRNA_matrix())
#     geneIds <- convert_to_HUGOIds(geneIds)
#     geneIds <- paste(geneIds, collapse=" ")
#     
#     #generate the HTML content
#     htmlContent <- sprintf(prefix, geneIds)
#     htmlContent
#   })
#   
#   #reactive value to store precomputed shiny results of mRNA data
#   mRNA_heatmap_compute_results <- reactiveValues() 
#   
#   mRNA_cache_time <- reactiveValues()
#   output$mRNA_cache_time = renderPrint({
#     print(mRNA_cache_time$time)
#   })
#   
#   output$microRNA_compute_time = renderPrint({
#     print(microRNA_heatmap_compute_results$time)
#   })  
