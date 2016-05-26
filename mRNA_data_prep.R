###
#get the AMP-AD data
###
# fetch without hgnc_symbol(inconsistencies in the MSBB data)
logFC = downloadFile('syn6041639') %>% dplyr::select(-hgnc_symbol)

# # connect to biomart and get hgnc from ensembl_gene_id
# ensembl <- useMart("ensembl")
# ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
# hg19_annot <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), 
#                     filters = "ensembl_gene_id", 
#                     values = unique(logFC$ensembl_gene_id), mart = ensembl) %>% 
#   mutate(hgnc_symbol=toupper(hgnc_symbol))
# save(hg19_annot, file='hg19_annot.RData')

load('hg19_annot.RData')

logFC <- left_join(logFC, hg19_annot, by="ensembl_gene_id")


# Filter by adjust p value < 1 and logFC greater than 
logFCFiltered <- logFC %>% filter(adj.P.Val <= 1, abs(logFC) >= 0)

#data frame of study, comparison used, ensemble_gene_ids, and logFC
#sorted by ensemble_gene_id
ad_data <-  logFCFiltered %>%
  dplyr::select(DataSetName, Comparison, ensembl_gene_id, logFC) %>%
  unite(DataSet.Comparison, DataSetName, Comparison, sep = ' ') %>%
  spread(DataSet.Comparison, logFC) %>%
  arrange(ensembl_gene_id)

# make ensembl_gene_ids rownames
ad_data_matrix <- ad_data %>% dplyr::select(-ensembl_gene_id) 
rownames(ad_data_matrix) <- ad_data$ensemble_gene_id

ad_data_matrix <- ad_data_matrix[,order(colnames(ad_data_matrix))]


phenoData <- logFCFiltered %>%
  dplyr::select(DataSetName, Comparison) %>%
  unite(DataSet.Comparison, DataSetName, Comparison, sep = ' ', remove = FALSE) %>% unique() %>%
  arrange(DataSet.Comparison)

rownames(phenoData) <- phenoData$DataSet.Comparison

featureData <- logFCFiltered %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol) %>% unique() %>%
  arrange(ensembl_gene_id)

rownames(featureData) <- featureData$ensembl_gene_id

eset.mRNA <- ExpressionSet(assayData=as.matrix(ad_data_matrix),
                           phenoData=AnnotatedDataFrame(phenoData), 
                           featureData=AnnotatedDataFrame(featureData))


# data frame for p-value
ad_data_pvalue <-  logFCFiltered %>%
  dplyr::select(DataSetName, Comparison, ensembl_gene_id, adj.P.Val, logFC) %>%
  unite(DataSet.Comparison, DataSetName, Comparison, sep = ' ', remove=FALSE)






