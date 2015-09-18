# get the PCBC samples raw miRNA counts
flog.info('Reading the PCBC raw miRNA Exp data from Synapse', name='synapse')

miRNA_normCounts <- synGet('syn2701942')
miRNA_normCounts <- read.delim(miRNA_normCounts@filePath, header=T, sep='\t', 
                               as.is=T, stringsAsFactors = F, check.names=F)

temp_rownames <- tolower(miRNA_normCounts$mir)
miRNA_normCounts$mir <- NULL
# miRNA_normCounts <- apply(miRNA_normCounts,2, function(x) as.numeric(x))
rownames(miRNA_normCounts) <- temp_rownames

# get the miRNA to genes mapping table from synapse
flog.info('Reading the miRNA to genes mapping table from Synapse', name='synapse')
miRNA_to_genes_obj <- synGet('syn3461627')
miRNA_to_genes <- fread(miRNA_to_genes_obj@filePath, header = FALSE, data.table=FALSE)
colnames(miRNA_to_genes) <- c("mirName", "ensembl_gene_id")

flog.info('Reading the miRNA metadata table from Synapse', name='synapse')
miRNAQuery <- sprintf("select %s from syn3219876",
                     paste(c(metadataIdCol, metadataColsToUse), collapse=","))
miRNAMetadataTable <- synTableQuery(miRNAQuery)
miRNA_metadata <- miRNAMetadataTable@values
miRNA_metadata <- unique(miRNA_metadata)
rownames(miRNA_metadata) <- miRNA_metadata[, metadataIdCol]
miRNA_metadata[, metadataIdCol] <- NULL

## Only keep samples in both
mirna_in_common <- intersect(rownames(miRNA_metadata), colnames(miRNA_normCounts))
miRNA_metadata <- miRNA_metadata[mirna_in_common, ]
miRNA_normCounts <- miRNA_normCounts[, mirna_in_common]

miRNA_features <- data.frame(explicit_rownames=rownames(miRNA_normCounts))
rownames(miRNA_features) <- rownames(miRNA_normCounts)

eset.miRNA <- ExpressionSet(assayData=as.matrix(miRNA_normCounts),
                            phenoData=AnnotatedDataFrame(miRNA_metadata),
                            featureData=AnnotatedDataFrame(miRNA_features))

## Add separate miRNA names to the featureData
foo <- fData(eset.miRNA) %>% 
  separate(explicit_rownames, c("miRNA1", "miRNA2"),
           sep=",", extra="error", remove=FALSE)

rownames(foo) <- rownames(fData(eset.miRNA))
fData(eset.miRNA) <- foo