#Load libraries
library(scater)
library(tidyverse)
library(Seurat)
library(biomaRt)

#Load full seurat data
parse <- readRDS("data/rds_objects/seurat_obj_before_QC_240201.rds")
multiome <- readRDS("data/rds_objects/meyer_scrna_seq_multiome_subset_with_celltypes.rds")


#Convert to sce
sce_parse <- as.SingleCellExperiment(parse)
sce_multiome <- as.SingleCellExperiment(multiome)

##Annotate by gene name

ensembl.parse <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", version = 99)
genemap_parse <- getBM(attributes=c('ensembl_gene_id', "external_gene_name", "gene_biotype"),
                 filters = 'external_gene_name',
                 values = rownames(sce_parse),
                 mart = ensembl.parse) %>%
  dplyr::rename(gene_id = ensembl_gene_id,
                symbol = external_gene_name) # %>%
#  distinct()
dim(genemap_parse)

ensembl.98 <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", version = 98)
genemap_multiome.98 <- getBM(attributes=c('ensembl_gene_id', "external_gene_name","gene_biotype"),
                       filters = 'external_gene_name',
                       values = rownames(sce_multiome),
                       mart = ensembl.98) %>%
  dplyr::rename(gene_id = ensembl_gene_id,
                symbol = external_gene_name) # %>%
#  distinct()
dim(genemap_multiome.98)

length(rownames(parse))
length(rownames(multiome))
length(genemap_multiome[,1])
length(genemap_parse[,1])
length(intersect(genemap_parse[,1], genemap_multiome[,1]))
length(intersect(genemap_parse[,2], genemap_multiome[,2]))
length(setdiff(genemap_parse[,1], genemap_multiome[,1]))
length(setdiff(genemap_multiome[,2], genemap_parse[,2]))

parse_unique <- setdiff(genemap_parse[,2], genemap_multiome[,2])
multiome_unique <- setdiff(genemap_multiome[,2], genemap_parse[,2])

#1 extract unique genes from each dataset
#2 check gene_biotype composition %




# featureData <- tibble(symbol=rownames(sce)) %>%
#   left_join(genemap, by="symbol") %>%
#   drop_na() %>%
#   filter(duplicated(symbol) == FALSE) %>%
#   filter(symbol != "")
# symbol <- dplyr::select(featureData, symbol) 
# gene_id <- dplyr::select(featureData, gene_id) 

sce_final <- 
  subset(sce,
         rownames(sce) %in% symbol$symbol,
  )

mcols(sce_final) <- 
  DataFrame(mcols(sce_final), featureData)
rownames(sce_final) <- rowData(sce_final)$gene_id

#convert back to seurat
seu_final <- as.Seurat(sce_final)


# ==================

genemap_parse |>
  filter(genemap_parse[,2] %in% setdiff(genemap_multiome[,2], genemap_parse[,2]))


