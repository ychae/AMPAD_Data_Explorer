library("dplyr")
library("org.Hs.eg.db")

#get the ensembl id's for the selected genes
convert_to_ensemblIds <- function(genes){
  filtered_df <- hg19_annot %>%
    filter(SYMBOL %in% genes  | ENSEMBL %in% genes |
             ENSEMBLTRANS %in% genes  | ENTREZID %in% genes )
  na.omit(unique(c(filtered_df$ENSEMBL)))
}

#get the HUGO ids for the selected genes
convert_to_HUGOIds <- function(genes){
  filtered_df <- hg19_annot %>%
    filter(SYMBOL %in% genes  | ENSEMBL %in% genes |
             ENSEMBLTRANS %in% genes  | ENTREZID %in% genes )
  na.omit(unique(c(filtered_df$SYMBOL)))
}


#get the ENTREZ ids for the selected genes
convert_to_EntrezIds <- function(genes){
  filtered_df <- hg19_annot %>%
    filter(SYMBOL %in% genes  | ENSEMBL %in% genes |
             ENSEMBLTRANS %in% genes  | ENTREZID %in% genes )
  na.omit(unique(c(filtered_df$ENTREZID)))
}
