library(DT)
library(ComplexHeatmap)

sample_gene_list <- c("APOE", "CD33", "CD44", "APP", "PER1", 
                      "PICALM", "VGF", "MAPT", "BIN1", "CD47", 
                      "TYROBP", "DOCK2", "FCER1G", "FYN", "MEF2C",
                      "NPTX2", "SNAP29", "STX4", "SYT1", "TARDBP")


#Define the server the logic
shinyServer(
  
  function(input, output, session) {
    session$sendCustomMessage(type="readCookie",
                              message=list(name='org.sagebionetworks.security.user.login.token'))
    
    foo <- observeEvent(input$cookie, {
      
      # Synapse login
      synapseLogin(sessionToken=input$cookie)
      
      #get gene expression data
      source("gene_exprs_data_prep.R")
      combined_metadata <- pData(eset.mRNA)
      # Sample column required for expression matrix filtering
      combined_metadata$Sample <- rownames(combined_metadata)
      
      dataset <- reactive({
        eset.mRNA
      })
      
      output$plotHelp <- renderUI({
        filter_type_text <- filter_type_help()
        p(class = "text-muted", filter_type_text)
      })
      
      output$metadataui <- renderUI({
        tagList(tags$table(class="table table-condensed",
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
                           
        ))
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
                                                       min=0.000001, max=0.05, value=0.00001, step=0.001),
                                           sliderInput('logFC', label=h6('Log Fold Change'),
                                                       min=0, max=4, value=2, step=0.1),
                                           actionButton("refreshValue", "Refresh")
                                           
                                           
               ))
      })
      
      filter_type_help <- reactive({
        "Plotting selected genes."
      })  
      
      filtered_dataset <- reactive({
        
        ds <- dataset()
        ds_filtered <- filter_by_metadata(input, ds)
        flog.debug(sprintf("filtered ds dims: %s", dim(ds_filtered)), name="server")
        
        user_feats <- user_submitted_features()
        
        if(input$custom_search == "Significant_Genes") {  
          adjPVal <- isolate(input$adjPVal)
          logFCThreshold <- isolate(input$logFC)
          
          flog.debug(sprintf('pval thresh = %s, logfc thresh = %s', adjPVal, logFCThreshold), name='server')
          filtered_by_pvalue <- rownames(as.matrix(exprs(eset.pval)))[rowSums(as.matrix(exprs(eset.pval)) <= adjPVal, na.rm = T) >= 1]
          filtered_by_logFC <- rownames(as.matrix(exprs(eset.logFC)))[rowSums(as.matrix(abs(exprs(eset.logFC))) >= log2(logFCThreshold), na.rm = T) >= 1]
          
          filtered_values <- intersect(filtered_by_pvalue, filtered_by_logFC)
          
          ds_filtered <- ds_filtered[unique(filtered_values)] #$ensembl_gene_id),]
          flog.debug(sprintf("filtered ds dims after pval: %s", dim(ds_filtered)), name="server")
        }
        
        else if ((input$custom_search == "Gene_List")) {
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
          hm <- heatmap_cache$heatmap_logfc
          
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
      
      heatmap_cache <- reactiveValues()
      
      #return the heatmap plot
      output$heatmap_gene <- renderPlot({  
        flog.debug("Making heatmap", name='server')
        m_eset <- filtered_dataset()
        m <- exprs(m_eset)
        
        rownames(m) <- fData(m_eset)$hgnc_symbol
        filtered_metadata <- pData(m_eset)
        print(length(rownames(m)))
        
        # Scale counts
        m = scale(m)
        m = t(scale(t(m)))
        
        ha <- HeatmapAnnotation(filtered_metadata[,-(1)], col = list(Study=c('MAYO'='palegreen2', 'MSBB'='yellow', 'ROSMAP'='violet'), 
                                                                     Status=c('AD'='aquamarine', 'Control'='chocolate', 'ND'='firebrick1', 'SD'='darkorchid', 'NCI'='thistle')) )
        
        if(length(rownames(m)) <= 20) {
          h <- memoise(Heatmap(m, top_annotation = ha, name = '', 
                               show_row_names = T, show_column_names = F, 
                               cluster_columns = F, show_column_dend = F))
        } else {
          h <- memoise(Heatmap(m, top_annotation = ha, name = '', 
                               show_row_names = F, show_column_names = F, 
                               cluster_columns = F, show_column_dend = F))  
        }
        heatmap_cache$heatmap_gene <- print(h)
        
      }
      )
      
      # show who is logged in
      output$loggedin <- renderUI({
        h6(paste("Logged in as:", synGetUserProfile()@displayName))
      })
    }
    )

})



