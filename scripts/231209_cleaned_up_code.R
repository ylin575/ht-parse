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
saveRDS(hthy, file = "data/rds_objects/seurat_obj_before_QC_240101.rds")


# Cell QC: add mitochondrial percentage to metadata

hthy[["percent.mt"]] <- PercentageFeatureSet(hthy, pattern = "^MT-")


# Cell QC: plots
plot <- VlnPlot(hthy, pt.size = 0,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, group.by = "sample")
plot
# SaveFigure(plot, "vln_QC", width = 12, height = 6)

# plot1 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# SaveFigure((plot1 + plot2),"scatter_QC", width = 12, height = 6, res = 200)

# feature scatter
plot1 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        group.by = "sample")
plot2 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        group.by = "sample")

plot1 + plot2
# Visualize QC metrics as histograms, find best cut off point

# nFeature_RNA
hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="Ht9"], 
     breaks=500, xlim=c(0,5000),
     xlab = "nFeature_RNA", ylab = "Frequency", main = "HT9")
hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="HT11"], 
     breaks=500, xlim=c(0,5000),
     xlab = "nFeature_RNA", ylab = "Frequency", main = "HT11")
hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="Ht12"], 
     breaks=500, xlim=c(0,5000),
     xlab = "nFeature_RNA", ylab = "Frequency", main = "HT12")
hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="Ht14"], 
     breaks=300, xlim=c(0,5000),
     xlab = "nFeature_RNA", ylab = "Frequency", main = "HT14")

# nCount_RNA
hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="Ht9"], 
     breaks=500, xlim=c(0,20000),
     xlab = "nCount_RNA", ylab = "Frequency", main = "HT9")
hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="HT11"], 
     breaks=500, xlim=c(0,20000),
     xlab = "nCount_RNA", ylab = "Frequency", main = "HT11")
hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="Ht12"], 
     breaks=500, xlim=c(0,20000),
     xlab = "nCount_RNA", ylab = "Frequency", main = "HT12")
hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="Ht14"], 
     breaks=500, xlim=c(0,20000),
     xlab = "nCount_RNA", ylab = "Frequency", main = "HT14")

# percent.mt
hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="Ht9"], 
     breaks=200, xlim=c(0,60),
     xlab = "percent.mt", ylab = "Frequency", main = "HT9")
hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="HT11"], 
     breaks=200, xlim=c(0,60),
     xlab = "percent.mt", ylab = "Frequency", main = "HT11")
hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="Ht12"], 
     breaks=200, xlim=c(0,60),
     xlab = "percent.mt", ylab = "Frequency", main = "HT12")
hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="Ht14"], 
     breaks=200, xlim=c(0,60),
     xlab = "percent.mt", ylab = "Frequency", main = "HT14")

# Compare QC metrics across samples using boxplot
boxplot(hthy@meta.data$nCount_RNA ~ hthy@meta.data$sample, 
          xlab = "sample #", ylab = "nCount_RNA")
boxplot(hthy@meta.data$nFeature_RNA ~ hthy@meta.data$sample, 
        xlab = "sample #", ylab = "nFeature_RNA")
boxplot(hthy@meta.data$percent.mt ~ hthy@meta.data$sample, 
        xlab = "sample #", ylab = "percent.mt")

#violin plot in linear scale
VlnPlot(hthy, pt.size = 0, group.by = "sample",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
#violin plot in log scale
VlnPlot(hthy, pt.size = 0, group.by = "sample", log = TRUE,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)


# check dimension of the original hthy object
dim(hthy)
dim(hthy@assays)
dim(hthy@meta.data)


# subset and calculate mean and 3*SD for percent.mt and nFeature_RNA, then
# downsample

ht9 <- subset(x = hthy, subset = sample == "Ht9")
mean.percent.mt <- mean(ht9@meta.data$percent.mt)
sd.percent.mt <- sd(ht9@meta.data$percent.mt)
upper.mt <- mean.percent.mt + 3 * sd.percent.mt
upper.feature <- mean(ht9@meta.data$nFeature_RNA) + 3*sd(ht9@meta.data$nFeature_RNA)
ht9_sub <- subset(ht9, subset = percent.mt <= upper.mt)
ht9_sub <- subset(ht9_sub, subset = nFeature_RNA >= 500 & nFeature_RNA <+ upper.feature)


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

# check original dimensions of the subsetted objects
dim(ht9)
dim(ht11)
dim(ht12)
dim(ht14)

dim(ht9_sub)
dim(ht11_sub)
dim(ht12_sub)
dim(ht14_sub)
# based on the dimensions above, determined that ht12 has 8032 cells, the lowest

# downsample to 8032 cells for each sample
ht9_sub <- subset(ht9_sub, downsample = 8032)
ht11_sub <- subset(ht11_sub, downsample = 8032)
ht12_sub <- subset(ht12_sub, downsample = 8032)
ht14_sub <- subset(ht14_sub, downsample = 8032)

# merge the filtered and downsampled samples back into one seurat object
hthy <- merge(x = ht9_sub, y = list(ht11_sub, ht12_sub, ht14_sub))

# check dimensions of hthy
dim(hthy)

# rejoin the layers after merging, this step seems to take a few minutes
hthy[["RNA"]] <- JoinLayers(hthy[["RNA"]])
# after rejoining layers, size of hthy is 913.1 MB

# save new seurat object into a rds object
saveRDS(hthy, file = "data/rds_objects/seurat_obj_after_QC_downsample_231209.rds")

# clean up workspace
rm(ht9, ht11, ht12, ht14, ht9_sub, ht11_sub, ht12_sub, ht14_sub, mat, mat_t)


# subset pbmc to remove low quality cells, multiplets, and high mt contaminations
# rna <- subset(rna, subset = nFeature_RNA > 500 & nFeature_RNA < nFeature_RNA_cutoff.upper & percent.mt < percent.mt.upper)





################################################




# re-run qc plots

# Cell QC: plots
plot <- VlnPlot(hthy, pt.size = 0,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, group.by = "sample")
plot
# SaveFigure(plot, "vln_QC", width = 12, height = 6)

# plot1 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# SaveFigure((plot1 + plot2),"scatter_QC", width = 12, height = 6, res = 200)

# feature scatter
plot1 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        group.by = "sample")
plot2 <- FeatureScatter(hthy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        group.by = "sample")

plot1 + plot2
# Visualize QC metrics as histograms, find best cut off point

# nFeature_RNA
hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="Ht9"], 
     breaks=500, xlim=c(0,5000),
     xlab = "nFeature_RNA", ylab = "Frequency", main = "HT9")
hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="HT11"], 
     breaks=500, xlim=c(0,5000),
     xlab = "nFeature_RNA", ylab = "Frequency", main = "HT11")
hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="Ht12"], 
     breaks=500, xlim=c(0,5000),
     xlab = "nFeature_RNA", ylab = "Frequency", main = "HT12")
hist(hthy@meta.data$nFeature_RNA[hthy@meta.data$sample=="Ht14"], 
     breaks=300, xlim=c(0,5000),
     xlab = "nFeature_RNA", ylab = "Frequency", main = "HT14")

# nCount_RNA
hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="Ht9"], 
     breaks=500, xlim=c(0,20000),
     xlab = "nCount_RNA", ylab = "Frequency", main = "HT9")
hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="HT11"], 
     breaks=500, xlim=c(0,20000),
     xlab = "nCount_RNA", ylab = "Frequency", main = "HT11")
hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="Ht12"], 
     breaks=500, xlim=c(0,20000),
     xlab = "nCount_RNA", ylab = "Frequency", main = "HT12")
hist(hthy@meta.data$nCount_RNA[hthy@meta.data$sample=="Ht14"], 
     breaks=500, xlim=c(0,20000),
     xlab = "nCount_RNA", ylab = "Frequency", main = "HT14")

# percent.mt
hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="Ht9"], 
     breaks=200, xlim=c(0,60),
     xlab = "percent.mt", ylab = "Frequency", main = "HT9")
hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="HT11"], 
     breaks=200, xlim=c(0,60),
     xlab = "percent.mt", ylab = "Frequency", main = "HT11")
hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="Ht12"], 
     breaks=200, xlim=c(0,60),
     xlab = "percent.mt", ylab = "Frequency", main = "HT12")
hist(hthy@meta.data$percent.mt[hthy@meta.data$sample=="Ht14"], 
     breaks=200, xlim=c(0,60),
     xlab = "percent.mt", ylab = "Frequency", main = "HT14")

# Compare QC metrics across samples using boxplot
boxplot(hthy@meta.data$nCount_RNA ~ hthy@meta.data$sample, 
        xlab = "sample #", ylab = "nCount_RNA")
boxplot(hthy@meta.data$nFeature_RNA ~ hthy@meta.data$sample, 
        xlab = "sample #", ylab = "nFeature_RNA")
boxplot(hthy@meta.data$percent.mt ~ hthy@meta.data$sample, 
        xlab = "sample #", ylab = "percent.mt")

#violin plot in linear scale
VlnPlot(hthy, pt.size = 0, group.by = "sample",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
#violin plot in log scale
VlnPlot(hthy, pt.size = 0, group.by = "sample", log = TRUE,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)




################################################




# normalize the data
hthy <- NormalizeData(hthy, normalization.method = "LogNormalize", scale.factor
                     = 10000)

# find highly variable features, these will be used for the PCA
# this function directly models the mean-variance relationship inherent in
# single cell data. Default returns 2000 genes
hthy <- FindVariableFeatures(hthy, selection.method = "vst", nfeatures = 6000)

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(hthy), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hthy)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2

# scaling using default # of features, i.e. from FindVariableFeatures
hthy <- ScaleData(hthy)

# perform scalings using all genes (use too much ram)
# all.genes <- rownames(rna)
# rna <- ScaleData(rna, features = all.genes)



# scaling with regressing out mt percent, using SCTransform
# rna <- SCTransform(rna, vars.to.regress = "percent.mt")

# ScaleData() does 2 things
#     1. shift the expression of each gene, so that the mean expression across cells
#        is 0
#     2.scales the expression of each gene, so that the variance across cells is 1
#         otherwise the heatmap is going to be overpowered by highly expressed
#         genes

# optional? use var.to.regress argument to remove unwanted sources of variation
# pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt")




# perform linear dimensional reduction - run PCA
#> what PCA does is explaining the variance in the data
#>    PC1 explains the most variance, followed by PC2, than PC3, ..., so on
hthy <- RunPCA(hthy, features = VariableFeatures(object = hthy))

# save seurat object after running pca
saveRDS(hthy, file = "data/rds_objects/seurat_obj_after_run_pca_231209.rds")

# examine the PCA results in 4 ways

# 1. simply print out the principle components
print(hthy[["pca"]], dims = 1:5, nfeatures = 5)

# 2. 
VizDimLoadings(hthy, dims = 1:2, reduction = "pca")

# 3. scatter plot
DimPlot(hthy, reduction = "pca")

# 4. heat map
DimHeatmap(hthy, dims = 1:15, cells = 500, balanced = TRUE)



# determine dimensionality of the dataset
# seurat clusters cells based on their PCA scores
# methods: heatmap, JackStraw plot, Elbow plot

# JackStraw Plot (can take long time for big data sets)
# rna <- JackStraw(rna, num.replicate = 100)
# rna <- ScoreJackStraw(rna, dims = 1:20)
# JackStrawPlot(rna, dims = 1:20)
# 
# # Elbow plot
ElbowPlot(hthy, ndims = 50)


# cluster the cells

# first construct a KNN graph based on the euclidean distance in PCA space, and
# refine the edge weights between any two cells based on the shared overlap in
# their local neighborhoods (Jaccard similarity), using the FindNeighbors() function
hthy <- FindNeighbors(hthy, dims = 1:20)

# then cluster the cells by applying modularity optimization techniques
# setting the resolution parameter between 0.4-1.2 typically returns good results
hthy <- FindClusters(hthy, resolution = 0.5)

# look at cluster IDs of the first 5 cells
head(Idents(hthy), 5)

# run non-linear dimensional reduction (UMAP/tSNE)
hthy <- RunUMAP(hthy, dims = 1:20)

# plot UMAP
# note that you can set `label = TRUE` or use the LabelClusters function to help
# label individual clusters
DimPlot(hthy, reduction = "umap", label = TRUE)
DimPlot(hthy, reduction = "umap", label = TRUE, group.by = "sample")
DimPlot(hthy, reduction = "umap", label = TRUE, split.by = "sample")

# You can save the object at this point so that it can easily be loaded back in 
# without having to rerun the computationally intensive steps performed above, or
# easily shared with collaborators.
# saveRDS seems not able to create new folders, so must create 'output' folder another way
# saveRDS(pbmc, file = "output/230913 analysis/pbmc_tutorial.rds")
saveRDS(hthy, file = "data/rds_objects/seurat_obj_after_run_umap_231209.rds")


# =======================

# ways to finding markers that define each cluster

# find all markers of cluster 2
# clusters2.markers <- FindMarkers(hthy, ident.1 = 2, min.pct = 0.25)
# head(clusters2.markers, n = 10)
# 
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(hthy, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)
# 
# # find markers for every cluster compared to all remaining cells, report only the positive
# # ones
# hthy.markers <- FindAllMarkers(hthy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# hthy.markers.top20 <- hthy.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 20, order_by = avg_log2FC)

# list of lineage-defining genes
mTEC <- c("EPCAM","AIRE","HLA-DRA","CCL19","CCL21")
cTEC <- c("EPCAM","PRSS16","PSMB11","LY75")
# Lymphoid <- c("PTPRC","CD3E","CD4","CD8A","IL7R","S100A4","CD19","MS4A1","GNLY")
# Myeloid <- c("FCER1A","CST3","ITGAX","IGTAM","CD14")

# list of other genes of interest
keratin <- "KRT"
keratin_num <- 1:30
keratin <- expand.grid(keratin,keratin_num)
keratin <- paste(keratin[,1], keratin[,2], sep = "")


# # plot violin plot of gene expression in each cluster
# VlnPlot(hthy, features = c(mTEC))
# VlnPlot(hthy, features = c(cTEC))
# VlnPlot(hthy, features = c(Lymphoid))
# VlnPlot(hthy, features = c(Myeloid))
# 
# # plot gene expression heatmap on top of UMAP clusters
# FeaturePlot(hthy, features = mTEC)
# FeaturePlot(hthy, features = cTEC)
# FeaturePlot(hthy, features = Lymphoid)
# FeaturePlot(hthy, features = Myeloid)

# plot gene expression heatmap on top of UMAP clusters split by sample
FeaturePlot(hthy, features = "AIRE", label = TRUE, split.by = "sample")
FeaturePlot(hthy, features = "HLA-DRA", label = TRUE, split.by = "sample")
FeaturePlot(hthy, features = "PRSS16", label = TRUE, split.by = "sample")
FeaturePlot(hthy, features = "LY75", label = TRUE, split.by = "sample")
FeaturePlot(hthy, features = "PTPRC", label = TRUE, split.by = "sample")
FeaturePlot(hthy, features = "TTN", label = TRUE, split.by = "sample")
FeaturePlot(hthy, features = "XCL1", label = TRUE)
FeaturePlot(hthy, features = keratin, label = TRUE)

# save object after looking at some features
saveRDS(hthy, file = "data/rds_objects/seurat_obj_after_run_umap_231209.rds")

# # you can plot raw counts as well
# VlnPlot(hthy, features = mTEC, slot = "counts", log = TRUE)
# 
# 
# 
# # DoHeatmap() generates an expression heatmap for given cells and features. In 
# # this case, we are plotting the top 20 markers (or all markers if less than 20) 
# # for each cluster.
# hthy.markers %>%
#   group_by(cluster) %>%
#   top_n(n = 5, wt = avg_log2FC) -> top5
# DoHeatmap(hthy, features = top5$gene) + NoLegend()


################################################################################


# assigning cell type identity to clusters

# original.cluster.ids <- c("0",
#                      "1",
#                      "2", 
#                      "3",
#                      "4", 
#                      "5", 
#                      "6",
#                      "7", 
#                      "8", 
#                      "9",
#                      "10",
#                      "11")


# YL01.cluster.ids <- c("pre-TECs",
#                      "mTECs_1",
#                      "cTECs_1", 
#                      "cTECs_2",
#                      "cTECs_3", 
#                      "mTECs_2 (Aire+)", 
#                      "mTECs_3",
#                      "thymocytes or nurse cells?", 
#                      "mTECs_4 (proliferating)", 
#                      "9 = ?",
#                      "endothelial cells",
#                      "fibroblasts")
# 
# new.cluster.ids <- c("0mTECs_2 (Aire+)",
#                      "1cTECs_1",
#                      "2x", 
#                      "3x",
#                      "4x", 
#                      "5cTECs_2", 
#                      "6x",
#                      "7x", 
#                      "8x", 
#                      "9x",
#                      "10x",
#                      "11x")
# 
# 
# names(new.cluster.ids) <- levels(hthy)
# hthy <- RenameIdents(hthy, new.cluster.ids)
# DimPlot(hthy, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# 
# 
# # save final object
# saveRDS(hthy, file = "data/rds_objects/231206_renamed-clusters.rds")
# 
# sessionInfo()



