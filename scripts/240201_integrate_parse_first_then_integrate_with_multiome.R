# Loading libraries and setting location paths
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)

# remove all objects (variables, functions, etc.) from the current working 
# environment or session.
rm(list = ls())

# data_path <- "/volume-general/analysis/data/pbmc/"
# fig_path <- "/volume-general/analysis/figures/pbmc/seurat/"


# Below we've included a few convenience functions for saving images and 
# reading/writing your Seurat object to disk. When you're working with larger 
# datasets, it's usually a good idea to save your progress after computationally 
# intensive steps so you can back track if you wish to do so.

# Convenience functions
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}


# Reading in data

# After reading in the data we'll perform basic filtering a on our expression 
# matrix to remove low quality cells and uninformative genes. The parameter 
# "min_genes" will keep cells that have at least 300 genes, and similarly, 
# "min_cells" will keep genes that are expressed in at least 5 cells. 
# Note: Seurat version 4.1 includes a convenience function to read Parse data from 
# the DGE folder. If you would like to use this function, please skip the code
# block below and see the section "Reading in data with Seurat >= 4.1"


DGE_folder <- "rawdata/all-sample/DGE_filtered/"

mat <- readMM(paste0(DGE_folder, "count_matrix.mtx"))

cell_meta <- read.delim(paste0(DGE_folder, "cell_metadata.csv"),
                        stringsAsFactor = FALSE, sep = ",")
genes <- read.delim(paste0(DGE_folder, "all_genes.csv"),
                    stringsAsFactor = FALSE, sep = ",")

cell_meta$bc_wells <- make.unique(cell_meta$bc_wells, sep = "_dup")
rownames(cell_meta) <- cell_meta$bc_wells
genes$gene_name <- make.unique(genes$gene_name, sep = "_dup")

# Setting column and rownames to expression matrix
colnames(mat) <- genes$gene_name
rownames(mat) <- rownames(cell_meta)
mat_t <- t(mat)

# Remove empty rownames, if they exist
mat_t <- mat_t[(rownames(mat_t) != ""),]

# Seurat version 5 or greater uses "min.features" instead of "min.genes"
hthy <- CreateSeuratObject(mat_t, min.features = 100, min.cells = 2, 
                           meta.data = cell_meta)

# remove mat and mat_t
rm(mat, mat_t, cell_meta, genes)

# When we create our Seurat object the plate well numbers (column names in the 
# expression matrix) from the experiment will automatically be assigned to the 
# cell identity slot. In other words, the program assumes this is how we want to 
# initially classify our cells. In general, we would like to avoid this behavior 
# so there isn't a bias towards a particular cell class when removing outliers.

# Setting our initial cell class to a single type, this will change after 
# clustering. 
hthy@meta.data$orig.ident <- factor(rep("hthy", nrow(hthy@meta.data)))
Idents(hthy) <- hthy@meta.data$orig.ident

# copied below 2 lines from parse tutorial, doesn't work (can't find functions)
# SaveObject(hthy, "seurat_obj_before_QC")
# hthy <- ReadObject("seurat_obj_before_QC")
# use saveRDS instead of SaveObject
saveRDS(hthy, file = "data/rds_objects/240201_parse_samples_before_QC.rds")
hthy <- readRDS(file = "data/rds_objects/240201_parse_samples_before_QC.rds")


# Cell QC: add mitochondrial percentage to metadata

hthy[["percent.mt"]] <- PercentageFeatureSet(hthy, pattern = "^MT-")


# # Cell QC: plots
# plot <- VlnPlot(hthy, pt.size = 0,
#                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
#                 ncol = 3, group.by = "sample")
# plot
# # SaveFigure(plot, "vln_QC", width = 12, height = 6)
# 
# # plot1 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # plot2 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # SaveFigure((plot1 + plot2),"scatter_QC", width = 12, height = 6, res = 200)
# 
# # feature scatter
# plot1 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "percent.mt",
#                         group.by = "sample")
# plot2 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
#                         group.by = "sample")
# 
# plot1 + plot2
# # Visualize QC metrics as histograms, find best cut off point
# 
# # nFeature_RNA
# hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="Ht9"], 
#      breaks=500, xlim=c(0,5000),
#      xlab = "nFeature_RNA", ylab = "Frequency", main = "HT9")
# hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="HT11"], 
#      breaks=500, xlim=c(0,5000),
#      xlab = "nFeature_RNA", ylab = "Frequency", main = "HT11")
# hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="Ht12"], 
#      breaks=500, xlim=c(0,5000),
#      xlab = "nFeature_RNA", ylab = "Frequency", main = "HT12")
# hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="Ht14"], 
#      breaks=300, xlim=c(0,5000),
#      xlab = "nFeature_RNA", ylab = "Frequency", main = "HT14")
# 
# # nCount_RNA
# hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="Ht9"], 
#      breaks=500, xlim=c(0,20000),
#      xlab = "nCount_RNA", ylab = "Frequency", main = "HT9")
# hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="HT11"], 
#      breaks=500, xlim=c(0,20000),
#      xlab = "nCount_RNA", ylab = "Frequency", main = "HT11")
# hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="Ht12"], 
#      breaks=500, xlim=c(0,20000),
#      xlab = "nCount_RNA", ylab = "Frequency", main = "HT12")
# hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="Ht14"], 
#      breaks=500, xlim=c(0,20000),
#      xlab = "nCount_RNA", ylab = "Frequency", main = "HT14")
# 
# # percent.mt
# hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="Ht9"], 
#      breaks=200, xlim=c(0,60),
#      xlab = "percent.mt", ylab = "Frequency", main = "HT9")
# hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="HT11"], 
#      breaks=200, xlim=c(0,60),
#      xlab = "percent.mt", ylab = "Frequency", main = "HT11")
# hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="Ht12"], 
#      breaks=200, xlim=c(0,60),
#      xlab = "percent.mt", ylab = "Frequency", main = "HT12")
# hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="Ht14"], 
#      breaks=200, xlim=c(0,60),
#      xlab = "percent.mt", ylab = "Frequency", main = "HT14")
# 
# # Compare QC metrics across samples using boxplot
# boxplot(hthy@meta.data$nCount_RNA ~ hthy@meta.data$sample, 
#         xlab = "sample #", ylab = "nCount_RNA")
# boxplot(hthy@meta.data$nFeature_RNA ~ hthy@meta.data$sample, 
#         xlab = "sample #", ylab = "nFeature_RNA")
# boxplot(hthy@meta.data$percent.mt ~ hthy@meta.data$sample, 
#         xlab = "sample #", ylab = "percent.mt")
# 
# #violin plot in linear scale
# VlnPlot(hthy, pt.size = 0, group.by = "sample",
#         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
#         ncol = 3)
# #violin plot in log scale
# VlnPlot(hthy, pt.size = 0, group.by = "sample", log = TRUE,
#         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
#         ncol = 3)
# 
# 
# # check dimension of the original hthy object
# dim(hthy)
# dim(hthy@assays)
# dim(hthy@meta.data)


# subset and calculate mean and 3*SD for percent.mt (ignore lower bound since
# they tend to be less than 0) and nFeature_RNA, then downsample

ht9 <- subset(x = hthy, subset = sample == "Ht9")
mean.percent.mt <- mean(ht9@meta.data$percent.mt)
sd.percent.mt <- sd(ht9@meta.data$percent.mt)
upper.mt <- mean.percent.mt + 3 * sd.percent.mt
upper.feature <- mean(ht9@meta.data$nFeature_RNA) + 3*sd(ht9@meta.data$nFeature_RNA)
ht9_sub <- subset(ht9, subset = percent.mt <= upper.mt)
ht9_sub <- subset(ht9_sub, subset = nFeature_RNA >= 500 & nFeature_RNA <= upper.feature)


ht11 <- subset(x = hthy, subset = sample == "HT11")
mean.percent.mt <- mean(ht11@meta.data$percent.mt)
sd.percent.mt <- sd(ht11@meta.data$percent.mt)
upper.mt <- mean.percent.mt + 3 * sd.percent.mt
upper.feature <- mean(ht11@meta.data$nFeature_RNA) + 3*sd(ht11@meta.data$nFeature_RNA)
ht11_sub <- subset(ht11, subset = percent.mt <= upper.mt)
ht11_sub <- subset(ht11_sub, subset = nFeature_RNA >= 500 & nFeature_RNA <= upper.feature)


ht12 <- subset(x = hthy, subset = sample == "Ht12")
mean.percent.mt <- mean(ht12@meta.data$percent.mt)
sd.percent.mt <- sd(ht12@meta.data$percent.mt)
upper.mt <- mean.percent.mt + 3 * sd.percent.mt
upper.feature <- mean(ht12@meta.data$nFeature_RNA) + 3*sd(ht12@meta.data$nFeature_RNA)
ht12_sub <- subset(ht12, subset = percent.mt <= upper.mt)
ht12_sub <- subset(ht12_sub, subset = nFeature_RNA >= 500 & nFeature_RNA <= upper.feature)


ht14 <- subset(x = hthy, subset = sample == "Ht14")
mean.percent.mt <- mean(ht14@meta.data$percent.mt)
sd.percent.mt <- sd(ht14@meta.data$percent.mt)
upper.mt <- mean.percent.mt + 3 * sd.percent.mt
upper.feature <- mean(ht14@meta.data$nFeature_RNA) + 3*sd(ht14@meta.data$nFeature_RNA)
ht14_sub <- subset(ht14, subset = percent.mt <= upper.mt)
ht14_sub <- subset(ht14_sub, subset = nFeature_RNA >= 500 & nFeature_RNA <= upper.feature)

# # check original dimensions of the subsetted objects
# dim(ht9)
# dim(ht11)
# dim(ht12)
# dim(ht14)
# dim(ht9_sub)
# dim(ht11_sub)
# dim(ht12_sub)
# dim(ht14_sub)
# based on the dimensions above, determined that ht12 has 8032 cells, the lowest

# 240122 no downsampling for our analysis needs
# downsample to 8032 cells for each sample
# ht9_sub <- subset(ht9_sub, downsample = 8032)
# ht11_sub <- subset(ht11_sub, downsample = 8032)
# ht12_sub <- subset(ht12_sub, downsample = 8032)
# ht14_sub <- subset(ht14_sub, downsample = 8032)

# remove original objects
rm(ht9, ht11, ht12, ht14)


# test integrate entire dataset (ht9, ht11, ht12, ht14)
# merge and join layers first for ht9,11,12,14, then merge with sarah's dataset
# merge the filtered and downsampled samples back into one seurat object
hthy2 <- merge(x = ht9_sub, y = list(ht11_sub, ht12_sub, ht14_sub))


hthy2 <- NormalizeData(hthy2)
hthy2 <- FindVariableFeatures(hthy2)
hthy2 <- ScaleData(hthy2)
hthy2 <- RunPCA(hthy2, features = VariableFeatures(object = hthy2))

hthy2 <- FindNeighbors(hthy2, dims = 1:30)
hthy2 <- FindClusters(hthy2, resolution = 0.5, 
                      cluster.name = "unintegrated_clusters")

hthy2 <- RunUMAP(hthy2, dims = 1:30, 
                 reduction.name = "umap.unintegrated")
DimPlot(hthy2, reduction = "umap.unintegrated", group.by = "sample")


# perform integration on hthy2: the 4 qc'd parse samples

hthy2wi <- IntegrateLayers(object = hthy2, method = CCAIntegration, 
                           orig.reduction = "pca", 
                           new.reduction = "integrated.cca",
                           verbose = FALSE) # 240102 took 15m55s: integrate complete parse and complete multiome
                                            # 240122 took 1hr: integrate qc>merge>integrated and multiome
                                            # 240123 took 1hr4min
                                            # 240129 took 55m 21s: integrate qc>merge>integrated parse

# rejoin the layers after merging, this step seems to take a few minutes
hthy2wi[["RNA"]] <- JoinLayers(hthy2wi[["RNA"]])
# after rejoining layers of the 4 parse samples, size = 3.9 GB

# remove objects
rm(ht9_sub, ht11_sub, ht12_sub, ht14_sub)

# save integrated parse dataset post-qc
saveRDS(hthy2wi, file = "data/rds_objects/240201_parse_samples_only_post-qc_merged_integrated.rds")
hthy2wi <- readRDS(file = "data/rds_objects/240201_parse_samples_only_post-qc_merged_integrated.rds")


hthy2wi <- FindNeighbors(hthy2wi, reduction = "integrated.cca", dims = 1:30)
hthy2wi <- FindClusters(hthy2wi, resolution = 0.5)
hthy2wi <- RunUMAP(hthy2wi, dims = 1:30, reduction = "integrated.cca")

#plot graphs before combining sample IDs
DimPlot(hthy2wi, reduction = "umap", group.by = "sample")
DimPlot(hthy2wi, reduction = "umap", split.by = "sample")

FeaturePlot(hthy2wi, reduction = "umap", features = "AIRE", split.by = "sample")
FeaturePlot(hthy2wi, reduction = "umap", features = "PRSS16", split.by = "sample")
FeaturePlot(hthy2wi, reduction = "umap", features = "PSMB11", split.by = "sample")
FeaturePlot(hthy2wi, reduction = "umap", features = "LY75", split.by = "sample")



# load in sarah's reference data
reference <- readRDS(file = "data/rds_objects/meyer_scrna_seq_multiome_subset_with_celltypes.rds")


all_samples_ni <- merge(x = hthy2wi, y = reference)


all_samples_ni <- NormalizeData(all_samples_ni)
all_samples_ni <- FindVariableFeatures(all_samples_ni)
all_samples_ni <- ScaleData(all_samples_ni)
all_samples_ni <- RunPCA(all_samples_ni, features = VariableFeatures(object = all_samples_ni))

all_samples_ni <- FindNeighbors(all_samples_ni, dims = 1:30)
all_samples_ni <- FindClusters(all_samples_ni, resolution = 0.5, 
                      cluster.name = "unintegrated_clusters")
all_samples_ni <- RunUMAP(all_samples_ni, dims = 1:30, 
                 reduction.name = "umap.unintegrated")

DimPlot(all_samples_ni, reduction = "umap.unintegrated", group.by = "sample")
DimPlot(all_samples_ni, reduction = "umap.unintegrated", group.by = "batch")

# combine sample IDs then re-do dimplots above
all_samples_ni@meta.data <- all_samples_ni@meta.data |>
  mutate(samples = ifelse(is.na(sample), batch, sample))
DimPlot(all_samples_ni, reduction = "umap.unintegrated", group.by = "samples")



# perform integration of parse and multiome dataset
all_samples_wi <- IntegrateLayers(object = all_samples_ni, method = CCAIntegration, 
                           orig.reduction = "pca", 
                           new.reduction = "integrated.cca",
                           verbose = FALSE)
# 240201 took 19m 20s to integrate

# rejoin the layers after merging
all_samples_wi[["RNA"]] <- JoinLayers(all_samples_wi[["RNA"]])

# save integrated object after integrating parse and multiome datasets
saveRDS(all_samples_wi, file = "data/rds_objects/240201_parse_multiome_integrated.rds")
all_samples_wi <- readRDS(file = "data/rds_objects/240201_parse_multiome_integrated.rds")


all_samples_wi <- FindNeighbors(all_samples_wi, reduction = "integrated.cca", dims = 1:30)
all_samples_wi <- FindClusters(all_samples_wi, resolution = 0.5)
all_samples_wi <- RunUMAP(all_samples_wi, dims = 1:30, reduction = "integrated.cca")

# plot graphs post-integration of parse and multiome data
DimPlot(all_samples_wi, group.by = "samples")
DimPlot(all_samples_wi, reduction = "umap", group.by = "samples")
DimPlot(all_samples_wi, reduction = "umap", split.by = "samples")



