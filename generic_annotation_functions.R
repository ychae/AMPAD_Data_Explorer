library("dplyr")
library("org.Hs.eg.db")

#gene annotation
cat('Preparing the hg19 annotation df.....')
k <- keys(org.Hs.eg.db,keytype="SYMBOL")
hg19_annot <- select(org.Hs.eg.db, keys=k, columns=c("GENENAME","ALIAS", "ENSEMBL", "ENSEMBLTRANS", "ENTREZID"), keytype="SYMBOL")
hg19_grpd <- hg19_annot %>%
  group_by(ENSEMBL) %>%
  summarise(ALIAS = paste(unique(ALIAS),collapse=", "),
            SYMBOL = paste(unique(SYMBOL),collapse=", "),
            GENENAME = paste(unique(GENENAME),collapse=", "),
            ENTREZID = paste(unique(ENTREZID),collapse=", ")
  )
hg19_grpd <- as.data.frame(hg19_grpd)
cat('Done \n\n')

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
