###
#get the PCBC samples geneExp normalized counts
###

flog.info('Reading the PCBC normalized mRNA Exp data from Synapse', name='synapse')
mRNA_NormCounts <- synGet('syn2701943')

#read in the file
mRNA_NormCounts <- read.delim(mRNA_NormCounts@filePath, header=T, sep='\t',
                              as.is=T, stringsAsFactors = F, check.names=F)

## remove version from ENSEMBL ID
rownames(mRNA_NormCounts) <- gsub('\\..*', '',mRNA_NormCounts$gene_id)
mRNA_NormCounts$symbol <- NULL
mRNA_NormCounts$gene_id <- NULL
mRNA_NormCounts$locus <- NULL
#apply(mRNA_NormCounts,2,class)
#mRNA_NormCounts <- as.data.frame(apply(mRNA_NormCounts,2,as.numeric))
#rownames(mRNA_NormCounts)

###
#get the metadata from synapse for PCBC geneExp samples
###
flog.info('Reading the PCBC mRNA metadata from Synapse', name='synapse')

mRNA_metadata <- synGet('syn2731147')

mRNA_metadata <- read.delim(mRNA_metadata@filePath, header=T, sep='\t',
                            as.is=T, stringsAsFactors = F, check.names=F)

rownames(mRNA_metadata) <- mRNA_metadata[,'Decorated Name']

#keep only that metadata for samples which we have expression data
mRNA_metadata <- mRNA_metadata[rownames(mRNA_metadata) %in% colnames(mRNA_NormCounts),]

