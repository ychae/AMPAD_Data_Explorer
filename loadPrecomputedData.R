## Load precomputed data
## Should check in precomputed_data/preCompute.R to make sure there are no
## needed updates

## Load in the hg19 annotations and grouped annotations
hg19_annot <- readRDS("precomputed_data/precomputed_hg19_annot.RDS")
hg19_grpd <- readRDS("precomputed_data/precomputed_hg19_grpd.RDS")

#sample gene list of the user input area
df <- read.table("precomputed_data/pre_selected_genelist.txt",sep="\t")
sample_gene_list <- as.character(unique(df$V5))

sample_miRNAs <- c("hsa-mir-627", "hsa-mir-34c", "hsa-let-7g",
                   "hsa-mir-19a", "hsa-mir-342")

#get the list siginificant genes from comparative analysis in synapse
flog.info('Reading the precomputed significant genelist')
sigGenes_lists <- readRDS("precomputed_data/precomputed_sigGenes_lists.rds")

#########
#read the precomputed enriched pathway list
########
df_precomputed_enrichedPathways_in_geneLists = readRDS("precomputed_data/precomputed_enrichedPathways_in_geneLists.rds")
df_precomputed_enrichedPathways_in_geneLists$pathways_with_pvalue =  paste(df_precomputed_enrichedPathways_in_geneLists$pathways,
                                                                           '#p.adj_',
                                                                           format.pval(df_precomputed_enrichedPathways_in_geneLists$p.adj,digits=2),
                                                                           sep='')
#creating a list of list 
precomputed_enrichedPathways_in_geneLists = split(df_precomputed_enrichedPathways_in_geneLists$pathways_with_pvalue,
                                                  df_precomputed_enrichedPathways_in_geneLists$significant_gene_list_name)


#HACK
#For each geneList add another PATHWAY TYPE "ALL" which indicates use all the pathways for the shiny SERVER/UI
# in this case genes in all the enriched pathways will be shown on the heatmap
precomputed_enrichedPathways_in_geneLists <- lapply(precomputed_enrichedPathways_in_geneLists,function(x) { x[length(x)+1] = 'ALL'; x})
