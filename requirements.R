# Requirements for running this Shiny app

library("devtools")
library("shinyIncubator")
library("synapseClient")
library("gdata")
library("shiny")
library("digest")
library("dplyr")
library("memoise")
library("org.Hs.eg.db")

#source the heatmap code
source_url("expression_heatmap.R")

#source generic heatmap functions
source_url("generic_annotation_functions.R")
