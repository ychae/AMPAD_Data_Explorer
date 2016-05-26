## It is assumed your working directory is where this file is

# Clear R console screen output
cat("\014")  

# Load required libraries
library(CovariateAnalysis)
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

synapseLogin()

# Get the latest commit of the code from git
thisFileName <- 'collateDiffExp.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/Brain_Reg_Net", ref="branch", refName='metaAnal')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0('code/R/',thisFileName))

# Get protein coding genes from biomaRt
ensembl=useMart("ensembl", host="www.ensembl.org")
ensemblHSapiens = useDataset("hsapiens_gene_ensembl",mart=ensembl)
humanProteinCodingGenes = getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), 
                                filters='biotype', 
                                values=c('protein_coding'), 
                                mart=ensemblHSapiens)

#### Get all covariates
covariates.id = c('ROSMAP' = 'syn6114444', 
                  'MSBB' = 'syn5573110',
                  'MAYO' = 'syn6088521')
covariates = lapply(covariates.id, downloadFile)

covariates$ROSMAP = covariates$ROSMAP %>%
  dplyr::select(Sampleid, msex, cogdx) %>%
  dplyr::rename(ID = Sampleid, Gender = msex, Status = cogdx) %>%
  dplyr::mutate(Gender = factor(Gender, labels = c('1' = 'male','2'='female')),
                Status = factor(Status, labels = c('1' = 'NCI', '2' = 'MCI', '3' = 'MCI.CI', 
                                                   '4' = 'AD', '5' = 'AD.CI', '6' = 'OD')),
                BrainRegion = 'DLPFC') %>%
  dplyr::filter(Status %in% c('NCI', 'AD'))

covariates$MSBB = covariates$MSBB %>%
  dplyr::select(SampleId, SEX, BrainRegion.Dx) %>%
  tidyr::separate(BrainRegion.Dx, c('BrainRegion', 'Dx'), sep = '\\.') %>%
  dplyr::rename(ID = SampleId, Gender = SEX, Status = Dx) %>%
  dplyr::mutate(Gender = factor(Gender, labels = c('M' = 'male','F'='female')),
                Status = factor(Status)) %>%
  dplyr::filter(Status %in% c('ND', 'SD'))

covariates$MAYO = covariates$MAYO %>%
  dplyr::select(ID, Gender, BrainRegion.Diagnosis) %>%
  tidyr::separate(BrainRegion.Diagnosis, c('BrainRegion', 'Dx'), sep = '\\.') %>%
  dplyr::rename(Status = Dx) %>%
  dplyr::mutate(Gender = factor(Gender, labels = c('M' = 'male','F'='female')),
                Status = factor(Status)) %>%
  dplyr::filter(Status %in% c('Control', 'AD'))

covariates = rbindlist(covariates, use.names = T, fill = T, idcol = 'Study')

#### Get raw expression values
logcpm.id = c('ROSMAP' = 'syn6114447', 'MSBB' = 'syn6117295', 'MAYO' = 'syn6121652')
logcpm = lapply(logcpm.id, downloadFile) %>%
  lapply(function(x){
    x = dplyr::filter(x, ensembl_gene_id %in% humanProteinCodingGenes$ensembl_gene_id) %>%
      dplyr::select(-hgnc_symbol)
  }) %>%
  join_all(type = 'full')

# Get all differential expresison results from synapse
diffexp.id = c('ROSMAP' = 'syn6114453',
               'MSBB' = 'syn5609009',
               'Mayo' = 'syn6088525')

diffexp = lapply(diffexp.id, downloadFile) %>%
  rbindlist(idcol = 'Study', fill=T, use.names=T) %>%
  filter(ensembl_gene_id %in% humanProteinCodingGenes$ensembl_gene_id,
         adj.P.Val <= 0.05, Study == 'ROSMAP')

# Filter logcpm based on differential expression
logcpm = filter(logcpm, ensembl_gene_id %in% diffexp$ensembl_gene_id)

counts = data.frame(logcpm[,-(1)])
rownames(counts) = logcpm$ensembl_gene_id
counts[is.na(counts)] = 0

covariates = data.frame(covariates)
rownames(covariates) = covariates$ID

counts = counts[,intersect(rownames(covariates), colnames(counts))]
covariates = covariates[intersect(rownames(covariates), colnames(counts)),]

# Plot individual heatmaps for every study
pl = dlply(covariates, .(Study), .fun = function(x, counts){
  browser()
  c = counts[,rownames(x)]
  ha = HeatmapAnnotation(x[,-(2)])
  Heatmap(c, top_annotation = ha)
}, list(counts))