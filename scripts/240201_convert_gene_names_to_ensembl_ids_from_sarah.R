library(scater)
library(tidyverse)
library(Seurat)
library(biomaRt)

#Load full seurat data
seu <- readRDS("your_seu_object")

#Convert to sce
sce <- as.SingleCellExperiment(seu)

##Annotate by gene name
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

genemap <- getBM(attributes=c('ensembl_gene_id', "external_gene_name"),
                 filters = 'external_gene_name',
                 values = rownames(sce),
                 mart = ensembl) %>%
  dplyr::rename(gene_id = ensembl_gene_id,
                symbol = external_gene_name) %>%
  distinct()

featureData <- tibble(symbol=rownames(sce)) %>%
  left_join(genemap, by="symbol") %>%
  drop_na() %>%
  filter(duplicated(symbol) == FALSE) %>%
  filter(symbol != "")
symbol <- dplyr::select(featureData, symbol) 
gene_id <- dplyr::select(featureData, gene_id) 

sce_final <- 
  subset(sce,
         rownames(sce) %in% symbol$symbol,
  )

mcols(sce_final) <- 
  DataFrame(mcols(sce_final), featureData)
rownames(sce_final) <- rowData(sce_final)$gene_id

#convert back to seurat
seu_final <- as.Seurat(sce_final)

