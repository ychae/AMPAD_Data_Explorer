
ad the memoised version of pheatmap

# Requirements for running this Shiny app

options(stringsAsFactors = F)
options(warn=-1)

library("devtools")
# library("shinyIncubator")
library("synapseClient")
# library("gdata")
library("shiny")
# library("digest")
library("dplyr")
library("tidyr")
library("memoise")
# library("org.Hs.eg.db")
library("futile.logger")
library(Biobase)
library(data.table)
library(CovariateAnalysis)

# Set up logging
flog.threshold(DEBUG, name='server')
flog.threshold(DEBUG, name='ui')
flog.threshold(DEBUG, name='global')
flog.threshold(INFO, name='synapse')

#login to synapse
synapseLogin()

flog.debug("Starting App", name="server")

#source the heatmap code
#source("expression_heatmap.R")

#source generic heatmap functions
#source("generic_annotation_functions.R")

#get the global functions
source("global_functions.R")

###############################
## Load precomputed data
###############################

# source("loadPrecomputedData.R")

###############################
## Load synapse data
###############################
## This may break!
## Cache the data thats actually used
# save(list=c("combined_metadata", "eset.mRNA", "eset.miRNA", "eset.meth", "meth_to_gene", "miRNA_to_genes", "pathways_list"),
#      file="cached_data.RData")
# f <- File("cached_data.RData", parentId="syn4108202")
# o <- synStore(f)



#get the mRNA expression data
#source("mRNA_data_prep.R")

#get gene expression data
source("gene_exprs_data_prep.R")
combined_metadata <- pData(eset.mRNA)

# Sample column required for expression matrix filtering
combined_metadata$Sample <- rownames(combined_metadata)



