
# Requirements for running this Shiny app

options(stringsAsFactors = F)
options(warn=-1)

library("devtools")
library("synapseClient")
library("shiny")
library("dplyr")
library("tidyr")
library("memoise")
library("futile.logger")
library(Biobase)
library(data.table)

# Set up logging
flog.threshold(DEBUG, name='server')
flog.threshold(DEBUG, name='ui')
flog.threshold(DEBUG, name='global')
flog.threshold(INFO, name='synapse')

#login to synapse
# synapseLogin()

flog.debug("Starting App", name="server")

#get the global functions
source("global_functions.R")

# Sample column required for expression matrix filtering
combined_metadata$Sample <- rownames(combined_metadata)



