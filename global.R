#load the memoised version of pheatmap

# Requirements for running this Shiny app

options(stringsAsFactors = F)
options(warn=-1)

library("devtools")
library("shinyIncubator")
library("synapseClient")
library("gdata")
library("shiny")
library("digest")
library("dplyr")
library("tidyr")
library("memoise")
library("org.Hs.eg.db")
library("futile.logger")
library(Biobase)
# Set up logging
flog.threshold(DEBUG, name='server')
flog.threshold(DEBUG, name='ui')
flog.threshold(DEBUG, name='global')
flog.threshold(INFO, name='synapse')

#login to synapse
synapseLogin()

#source the heatmap code
source("expression_heatmap.R")

#source generic heatmap functions
source("generic_annotation_functions.R")

#get the global functions
source("global_functions.R")

###############################
## Load precomputed data
###############################

source("loadPrecomputedData.R")

###############################
## Load synapse data
###############################


## This may break!
## Cache the data thats actually used
# save(list=c("combined_metadata", "eset.mRNA", "eset.miRNA", "eset.meth", "meth_to_gene", "miRNA_to_genes", "pathways_list"),
#      file="cached_data.RData")
# f <- File("cached_data.RData", parentId="syn4108202")
# o <- synStore(f)

use_cache <- FALSE
cacheId <- "syn4108151"

if (use_cache) {
  ## Caching for testing
  o <- synGet(cacheId)
  load(getFileLocation(o))
  flog.debug("Using cached data loaded from Synapse", name="server")
} else {
  # Use only these metadata columns
  metadataColsToUse <- c("Cell_Line_Type", "Reprogramming_Gene_Combination", 
                         "Reprogramming_Vector_Type", "Tissue_of_Origin", "Diffname_short",
                         "Cell_Type_of_Origin", "Gender", "Originating_Lab_ID",
                         "Cell_Line_of_Origin", "Donor_ID", "Originating_Lab", "Cell_Type",
                         "Culture_Conditions")
  # metadataColsToUse <- c("Cell_Line_Type")
  metadataIdCol <- "UID"
  
  #get the MSigDB data
  source("msigdb_data_prep.R")
  
  #get the mRNA expression data
  source("mRNA_data_prep.R")
  
  #get the miRNA expression data
  source("miRNA_data_prep.R")
  
  #get the methylation data
  source("methylation_data_prep.R")
  
  #prepare single global metadata
  combined_metadata <- rbind(mRNA_metadata, miRNA_metadata, meth_metadata, deparse.level = 0)
  
  # Sample column required for expression matrix filtering
  combined_metadata$Sample <- rownames(combined_metadata)
}