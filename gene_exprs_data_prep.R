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

colnames(logcpm)[1] <- "GeneID"
# 
counts <- logcpm %>% 
  filter(!is.na(GeneID)) %>% 
#   dplyr::select(GeneID, ensembl_gene_id) %>%
#   tidyr::unite(GeneID.ensembl_gene_id, GeneID, ensembl_gene_id, sep = "_") %>% 
  data.frame
  #   spread(GeneID.ensembl_gene_id) %>%
#   arrange(ensembl_gene_id)


rownames(counts) <- counts$ensembl_gene_id
counts$ensembl_gene_id <- NULL
counts[is.na(counts)] <- 0

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


