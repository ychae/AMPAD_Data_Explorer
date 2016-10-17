## It is assumed your working directory is where this file is

# Clear R console screen output
cat("\014")  

# Load required libraries
#library(CovariateAnalysis)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

library(synapseClient)
library(knitr)
library(githubr)

library(biomaRt)
library(ComplexHeatmap)

library(shiny)
library(grid)
library(Biobase)
library(ggplot2)

# synapseLogin()

# Utility function to download tsv or csv file from synapse and load it in to memory
downloadFile <- function(id, ...){
  tmp = data.table::fread(synapseClient::synGet(id)@filePath, data.table=F, header=T, ...)
}

# Covariates
covariates <- downloadFile('syn6132532')

# logcpm 
logcpm <- downloadFile('syn6132534')

# diffexp
diffexp <- downloadFile('syn6132536')

# Filter logcpm based on differential expression
# counts <- logcpm %>%
#   dplyr::select(-Gene.ID) %>%
#   filter(ensembl_gene_id %in% diffexp$ensembl_gene_id) %>%
#   group_by(ensembl_gene_id) %>% 
#   slice(1) %>% 
#   ungroup() %>%
#   as.data.frame()

##############
# Heatmap data
##############

# Remove period from Gene.ID colname
colnames(logcpm)[1] <- "GeneID"

# Remove NA gene Ids
counts <- logcpm %>% 
  filter(!is.na(GeneID)) %>% 
  filter(ensembl_gene_id %in% diffexp$ensembl_gene_id) %>%
  as.data.frame

# Remove non-ensemble_gene_ids (e.g. '_alignment_not_sufficient_')
counts <- counts[grepl("ENSG", counts$ensembl_gene_id),]

# name the rows by ensembl gene id
rownames(counts) <- counts$ensembl_gene_id
# counts$ensembl_gene_id <- NULL
counts[is.na(counts)] <- 0

# df of ensemble gene ids to filter by
counts_genes <- data.frame(counts$ensembl_gene_id)
colnames(counts_genes) <- "ensembl_gene_id"


# Add rownames to covariates
covariates <- data.frame(covariates)
rownames(covariates) <- covariates$ID

# Filter and arrange covariates and counts
counts <- counts[,intersect(rownames(covariates), colnames(counts))]
covariates <- covariates[intersect(rownames(covariates), colnames(counts)),]

# Arrange covariates
covariates <- covariates %>%
  arrange(Study, BrainRegion, Status, Gender) %>%
  dplyr::select(ID, Study, BrainRegion, Status, Gender) %>%
  data.frame
rownames(covariates) = covariates$ID

# Convert factor variables to character 
i <- sapply(covariates, is.factor)
covariates[i] <- lapply(covariates[i], as.character)

#data frame of study, comparison used, ensemble_gene_ids, and logFC
#sorted by ensemble_gene_id
ad_data <-  diffexp %>%
  filter(!is.na(hgnc_symbol)) %>%
  filter(ensembl_gene_id %in% counts_genes$ensembl_gene_id) %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol, Study, Tissue, Contrast, logFC) %>%
  tidyr::unite(Study.Tissue.Contrast, Study, Tissue, Contrast, sep = '_') %>%
  spread(Study.Tissue.Contrast, logFC) %>%
  arrange(ensembl_gene_id)

# make ensembl_gene_ids rownames
ad_data_matrix <- ad_data %>% 
  dplyr::select(-ensembl_gene_id, -hgnc_symbol) 
rownames(ad_data_matrix) <- ad_data$ensembl_gene_id

phenoData <- diffexp %>%
  filter(!is.na(hgnc_symbol)) %>%
  dplyr::select(Study, Tissue, Contrast) %>%
  unite(Study.Tissue.Contrast, Study, Tissue, Contrast, sep = '_', remove = FALSE) %>% 
  unique() %>%
  arrange(Study.Tissue.Contrast) %>%
  data.frame

rownames(phenoData) <- phenoData$Study.Tissue.Contrast

featureData <- diffexp %>%
  filter(!is.na(hgnc_symbol)) %>%
  filter(ensembl_gene_id %in% counts_genes$ensembl_gene_id) %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol) %>%
  unique() %>%
  arrange(ensembl_gene_id) %>%
  data.frame

rownames(featureData) <- featureData$ensembl_gene_id

eset.logFC <- ExpressionSet(assayData=as.matrix(ad_data_matrix),
                            phenoData=AnnotatedDataFrame(phenoData), 
                            featureData=AnnotatedDataFrame(featureData))
# data frame for p-value
ad_data_pvalue <-  diffexp %>%
  filter(!is.na(hgnc_symbol)) %>%
  filter(ensembl_gene_id %in% counts_genes$ensembl_gene_id) %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol, Study, Tissue, Contrast, adj.P.Val) %>%
  tidyr::unite(Study.Tissue.Contrast, Study, Tissue, Contrast, sep = '_') %>%
  spread(Study.Tissue.Contrast, adj.P.Val) %>%
  arrange(ensembl_gene_id)

ad_data_matrix <- ad_data_pvalue %>% 
  dplyr::select(-ensembl_gene_id, -hgnc_symbol) 
rownames(ad_data_matrix) <- ad_data_pvalue$ensembl_gene_id

eset.pval <- ExpressionSet(assayData=as.matrix(ad_data_matrix),
                           phenoData=AnnotatedDataFrame(phenoData), 
                           featureData=AnnotatedDataFrame(featureData))

# Final output for heatmap viewer/explorer
eset.mRNA <- ExpressionSet(assayData=as.matrix(counts)[rownames(featureData), rownames(covariates)],
                           phenoData=AnnotatedDataFrame(covariates), 
                           featureData=AnnotatedDataFrame(featureData))



####################
# Split violin plots, effect size plots, and data table 
####################

# Filter logcpm based on differential expression
counts_vp <- logcpm %>%
  dplyr::select(-Gene.ID) %>%
  filter(ensembl_gene_id %in% diffexp$ensembl_gene_id) 

# Remove non-ensembl gene id values
counts_vp <- counts_vp[grepl("ENSG", counts_vp$ensembl_gene_id),]


library(reshape2)
logcpm_long <- melt(counts_vp,
                    id.vars = c("ensembl_gene_id"),
                    variable.name = "ID",
                    value.name = "expression_level")

covar_edited <- covariates %>%
  dplyr::select(-Gender)

# Merge logcpm with covariates by sample_id
covar_logcpm <- merge(x = logcpm_long, y = covar_edited, by = "ID", all.x = T)

# change column name of BrainRegion to Tissue to match colnames in diffexp
colnames(covar_logcpm)[6] <- "Tissue"

# Get the columns from diffexp that are needed
diffexp_data <-  diffexp %>%
  filter(!is.na(hgnc_symbol)) %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol, Study, Tissue, Contrast, logFC, adj.P.Val) %>%
  arrange(ensembl_gene_id)

# Merge diffexp_data and covar_logcpm
final_merged_df <- merge(x = covar_logcpm, y = diffexp_data, by = c("ensembl_gene_id", "Study", "Tissue"))

# Convert the variable Status to a factor variable
final_merged_df$Status <- as.factor(final_merged_df$Status)

# convert to numeric
final_merged_df$expression_level <- as.numeric(final_merged_df$expression_level)

# filter out any NA expression levels
final_data_df <- final_merged_df %>%
  dplyr::filter(!is.na(expression_level)) 

# pval and logfc as 3 sig figs
final_data_df$adj.P.Val <- lapply(final_data_df$adj.P.Val, signif, 3)
final_data_df$logFC <- sapply(final_data_df$logFC, signif, 3)

# subset of mayo data
mayo_df <- final_data_df[final_data_df$Study == "MAYO", ]

# subset of msbb data
msbb_df <- final_data_df[final_data_df$Study == "MSBB", ]

# subset of rosmap data
rosmap_df <- final_data_df[final_data_df$Study == "ROSMAP", ]

# Legends for abbreviations
brainRegion <- "*Brain Tissue Legend: CER=Cerebellum, TCX=Temporal Cortex, 
FP=Frontal Pole, IFG=Inferior Frontal Gyrus, PHG=Parahippocampal Gyrus, 
STG=Superior Temporal Gyrus, DLPFC=Dorsolateral Prefrontal Cortex"
diseaseStatus <- "+Disease Status Legend: AD=Alzheimer's Disease, ND=No Dementia, SD=Severe Dementia,
NCI=No Cognitive Impairment"


# Data table for data table column
data_table_tab <- diffexp_data %>%
  dplyr::select(-ensembl_gene_id) 

#Rename pval column to not have periods
colnames(data_table_tab)[6] <- "fdr_adj_p_val" 
#Space holder for when there are wiki pages for each
data_table_tab$synId <- NA 

#####################
# Split violin plots
#####################

split_violin_fx <- function(dataframe, gene_name) {
  gene_df <- filter(dataframe, hgnc_symbol == gene_name)
  
  labelfun <- function(variable, value) {
    p <- lapply(value, function(v) 
      as.character(unique(gene_df$adj.P.Val[which(gene_df[,variable] == v)]))
    )
    unlist(p)
    l <- lapply(value, function(v) 
      as.character(unique(gene_df$logFC[which(gene_df[,variable] == v)]))
    )
    unlist(l)
    
    return(paste(value,"*", "\n" ,"FDR: ", p, "\n", "Log FC: ", l))
  }
  
  pdat <- gene_df %>%
    mutate(Tissue=factor(Tissue)) %>%
    group_by(Tissue, Status) %>%
    do(data.frame(loc = density(.$expression_level)$x,
                  dens = density(.$expression_level)$y))
  
  pdat$dens2 <- ifelse(pdat$Status == 'Control' | pdat$Status == 'ND' | pdat$Status == 'NCI', pdat$dens * -1, pdat$dens)
  
  split_violin_plot <- ggplot(pdat, aes(dens2, loc, fill = Status, group = interaction(Status, Tissue))) + 
    geom_polygon() +
    facet_wrap(~ Tissue, nrow = 1, labeller = labelfun) +
    labs(title=gene_df$Study) + 
    theme_bw() +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.title.y=element_blank()) + 
    scale_fill_discrete(name = bquote("Status"^"+")) +
    xlim(-1.5,1.5)
  
  return(split_violin_plot)
}


split_violin_fx_1 <- function(dataframe, gene_name) {
  gene_df <- filter(dataframe, hgnc_symbol == gene_name)
  
  labelfun <- function(variable, value) {
    p <- lapply(value, function(v) 
      as.character(unique(gene_df$adj.P.Val[which(gene_df[,variable] == v)]))
    )
    unlist(p)
    l <- lapply(value, function(v) 
      as.character(unique(gene_df$logFC[which(gene_df[,variable] == v)]))
    )
    unlist(l)
    
    return(paste(value,"*", "\n" ,"FDR: ", p, "\n", "Log FC: ", l))
  }
  
  pdat <- gene_df %>%
    mutate(Tissue=factor(Tissue)) %>%
    group_by(Tissue, Status) %>%
    do(data.frame(loc = density(.$expression_level)$x,
                  dens = density(.$expression_level)$y))
  
  pdat$dens2 <- ifelse(pdat$Status == 'Control' | pdat$Status == 'ND' | pdat$Status == 'NCI', pdat$dens * -1, pdat$dens)
  
  split_violin_plot <- ggplot(pdat, aes(dens2, loc, fill = Status, group = interaction(Status, Tissue))) + 
    geom_polygon() +
    facet_wrap(~ Tissue, nrow = 1, labeller = labelfun) +
    labs(list(title=gene_df$Study, y="Expression Level")) + 
    theme_bw() +
    theme(legend.position = "bottom",
          axis.title.x = element_blank()) + 
    scale_fill_discrete(name = bquote("Status"^"+")) +
    xlim(-1.5,1.5) 
  
  return(split_violin_plot)
}


#################
###Temp Effect Size Plots Columns
#################

numRows <- dim(diffexp)[1]

randomNums <- runif(numRows)

diffexp_fp <- diffexp %>%
  mutate(upperLogfc = logFC + randomNums, lowerLogfc = logFC - randomNums)

fp_mayo <- diffexp_fp[diffexp_fp$Study == 'MAYO', ]
fp_msbb <- diffexp_fp[diffexp_fp$Study == 'MSBB', ]
fp_rosmap <- diffexp_fp[diffexp_fp$Study == 'ROSMAP', ]

##########
# Effect size plots
##########

# Factor tissue so plot is grouped together by Study
diffexp_fp$Tissue<-factor(diffexp_fp$Tissue, levels=c('DLPFC', 'STG', 'PHG', 'IFG', 'FP', 'TCX', 'CER'))


forest_plot_fx <- function(dataframe, gene_name) {
  
  gene_df <- filter(dataframe, hgnc_symbol == gene_name)
  
  forest_plot <- ggplot(gene_df, aes(x=Tissue, y=logFC, ymin=lowerLogfc, ymax=upperLogfc)) +
    theme_bw() +
    geom_pointrange(aes(color=Study)) +
    geom_hline(yintercept = 0, linetype = 2) +
    coord_flip() +
    xlab('Tissue') + 
    labs(title=gene_name) + 
    theme(legend.position='bottom')  
  return(forest_plot)
}



