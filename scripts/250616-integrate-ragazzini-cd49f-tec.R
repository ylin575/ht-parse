# goal: re-assess integration of ragazzini dataset with our own data

# Loading libraries and setting location paths
library(Seurat)
library(dplyr)
library(Matrix)

# load stem tec data
ctec.cd49f.data <- Read10X(data.dir = "stemcell/ctec_cd49f")
mtec.cd49f.data <- Read10X(data.dir = "stemcell/mtec_cd49f")
ctec.data <- Read10X(data.dir = "stemcell/ctec")
mtec.data <- Read10X(data.dir = "stemcell/mtec")

# create seurat objects
ctec.cd49f <- CreateSeuratObject(counts = ctec.cd49f.data, 
                                 project = "ctec.cd49f", 
                                 min.cells = 2, 
                                 min.features = 100)
mtec.cd49f <- CreateSeuratObject(counts = mtec.cd49f.data, 
                                 project = "mtec.cd49f", 
                                 min.cells = 2, 
                                 min.features = 100)
ctec <- CreateSeuratObject(counts = ctec.data, 
                           project = "ctec", 
                           min.cells = 2, 
                           min.features = 100)
mtec <- CreateSeuratObject(counts = mtec.data, 
                           project = "mtec", 
                           min.cells = 2, 
                           min.features = 100)

# perform QC

# calculate percent mitochondrial genes and add to the meta data
mt_pct <- function(x) {
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
}

mt_pct(ctec).  # doesn't work

ctec[["percent.mt"]] <- PercentageFeatureSet(ctec, pattern = "^MT-")
ctec.cd49f[["percent.mt"]] <- PercentageFeatureSet(ctec.cd49f, pattern = "^MT-")
mtec[["percent.mt"]] <- PercentageFeatureSet(mtec, pattern = "^MT-")
mtec.cd49f[["percent.mt"]] <- PercentageFeatureSet(mtec.cd49f, pattern = "^MT-")



# RNA feature, count, percent mitochondrial

# violin plots
vln_qc <- function(x) {
  plot1 <- VlnPlot(x, pt.size = 0,
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3)
  plot1
  # SaveFigure(plot, "vln_QC", width = 12, height = 6)
}

vln_qc(ctec)
vln_qc(ctec.cd49f)
vln_qc(mtec)
vln_qc(mtec.cd49f)


# feature scatter plots
feature_qc <- function(x) {
  plot1 <- FeatureScatter(x, 
                          feature1 = "nCount_RNA", 
                          feature2 = "percent.mt")
  plot2 <- FeatureScatter(x, 
                          feature1 = "nCount_RNA", 
                          feature2 = "nFeature_RNA")
  plot1 + plot2
}

feature_qc(ctec)
feature_qc(ctec.cd49f)
feature_qc(mtec)
feature_qc(mtec.cd49f)



# Visualize QC metrics as histograms, find best cut off point

qc_3 <- function(x) {
  # nFeature_RNA
  hist(x@meta.data$nFeature_RNA,
       breaks=500, xlim=c(0,10000),
       xlab = "nFeature_RNA", ylab = "Frequency", 
       main = x@meta.data$orig.ident[1])
  # nCount_RNA
  hist(x@meta.data$nCount_RNA,
       breaks=500, xlim=c(0,25000),
       xlab = "nCount_RNA", ylab = "Frequency", 
       main = x@meta.data$orig.ident[1])
  # percent.mt
  hist(x@meta.data$percent.mt,
       breaks=200, xlim=c(0,50),
       xlab = "percent.mt", ylab = "Frequency", 
       main = x@meta.data$orig.ident[1])
}

qc_3(ctec)
qc_3(mtec)
qc_3(ctec.cd49f)
qc_3(mtec.cd49f)


# check dimension of the original object
dim(ctec)
dim(ctec@assays)
dim(ctec@meta.data)


# write a function to qc filter out bad quality cells

# mean.percent.mt <- mean(ctec@meta.data$percent.mt)
# sd.percent.mt <- sd(ctec@meta.data$percent.mt)
# upper.mt <- mean.percent.mt + 3 * sd.percent.mt
# upper.feature <- mean(ctec@meta.data$nFeature_RNA) + 3*sd(ctec@meta.data$nFeature_RNA)
# ctec <- subset(ctec, subset = percent.mt <= upper.mt)
# ctec <- subset(ctec, subset = ctec@meta.data$nFeature_RNA >= 500 & ctec@meta.data$nFeature_RNA <= upper.feature)

filter <- function(x) {
  mean.percent.mt <- mean(x@meta.data$percent.mt)
  sd.percent.mt <- sd(x@meta.data$percent.mt)
  upper.mt <- mean.percent.mt + 3 * sd.percent.mt
  upper.feature <- mean(x@meta.data$nFeature_RNA) + 3*sd(x@meta.data$nFeature_RNA)
  x <- subset(x, subset = percent.mt <= upper.mt)
  x <- subset(x, subset = x@meta.data$nFeature_RNA >= 500 & x@meta.data$nFeature_RNA <= upper.feature)
}

ctec <- filter(ctec)
ctec.cd49f <- filter(ctec.cd49f)
mtec <- filter(mtec)
mtec.cd49f <- filter(mtec.cd49f)



# load in parse data with no cell type labels, do QC
parse <- readRDS(file = "data/rds_objects/240205_parse_dataset_before_QC.rds")
parse[["percent.mt"]] <- PercentageFeatureSet(parse, pattern = "^MT-")

ht9 <- subset(x = parse, subset = sample == "Ht9")
ht11 <- subset(x = parse, subset = sample == "HT11")
ht12 <- subset(x = parse, subset = sample == "Ht12")
ht14 <- subset(x = parse, subset = sample == "Ht14")

ht9 <- filter(ht9)
ht11 <- filter(ht11)
ht12 <- filter(ht12)
ht14 <- filter(ht14)

ht9@meta.data$orig.ident <- "HT9"
ht11@meta.data$orig.ident <- "HT11"
ht12@meta.data$orig.ident <- "HT12"
ht14@meta.data$orig.ident <- "HT14"



# laod in sarah's multioime data with cell type labels
reference <- readRDS(file = "data/rds_objects/meyer_snrna_seq_multiome_full_data_with_celltypes.rds")

# get the batch (sample) names from reference data set
batches <- unique(reference@meta.data$batch)

ref1 <- subset(x = reference, subset = batch == batches[1])
ref2 <- subset(x = reference, subset = batch == batches[2])
ref3 <- subset(x = reference, subset = batch == batches[3])

# rename orig.ident to batch id
ref1@meta.data$orig.ident <- ref1@meta.data$batch[1]
ref2@meta.data$orig.ident <- ref2@meta.data$batch[1]
ref3@meta.data$orig.ident <- ref3@meta.data$batch[1]

# merge all samples from parse, multiome, and the ragazzini cd49f tecs
hthy <- merge(x = ctec.cd49f, y = list(ht9, ht11, ht12, ht14, 
                                        ref1, ref2, ref3,
                                        mtec.cd49f))

saveRDS(hthy, file = "data/rds_objects/250615-post-merge.rds")


# standard analysis steps.     note: piping doesn't work. why?
hthy <- NormalizeData(hthy)
hthy <- FindVariableFeatures(hthy)
hthy <- ScaleData(hthy)
hthy <- RunPCA(hthy, features = VariableFeatures(object = hthy))

# check elbow plot
ElbowPlot(hthy)

hthy <- FindNeighbors(hthy, dims = 1:15)
hthy <- FindClusters(hthy, resolution = 0.6, 
                     cluster.name = "unintegrated_clusters")
hthy <- RunUMAP(hthy, dims = 1:15, 
                reduction.name = "umap.unintegrated")

# look at unintegrated umap plot by orig.ident
DimPlot(hthy, reduction = "umap.unintegrated", group.by = "orig.ident")
DimPlot(hthy, reduction = "umap.unintegrated", split.by = "orig.ident")


#save before integration
saveRDS(hthy, file = "data/rds_objects/250615_before_integration.rds")
#hthy <- readRDS("data/rds_objects/250615_before_integration.rds")


# perform integration on the merged samples; wi = with integration
hthy <- IntegrateLayers(object = hthy, method = CCAIntegration, 
                               orig.reduction = "pca", 
                               new.reduction = "integrated.cca",
                               verbose = TRUE) # 250615 took 1h 43m 35s

# save object after integrating, before joining layers
saveRDS(hthy, file = "data/rds_objects/250615-after-integration.rds")
hthy <- readRDS("data/rds_objects/250615-after-integration.rds")


# rejoin the layers after integration
hthy[["RNA"]] <- JoinLayers(hthy[["RNA"]])

# analyze post integration
hthy <- FindNeighbors(hthy, reduction = "integrated.cca", dims = 1:15)
hthy <- FindClusters(hthy, resolution = 0.6)
hthy <- RunUMAP(hthy, dims = 1:15, reduction = "integrated.cca",
                n.components = 3L)  # added n.components for plotting 3D UMAP
                                    # see code below

# combine sample IDs
#hthy2wi@meta.data <- hthy2wi@meta.data |>
#  mutate(samples = ifelse(is.na(sample), batch, sample))

# look at integrated umap plot by orig.ident
DimPlot(hthy, reduction = "integrated.cca", group.by = "orig.ident")
DimPlot(hthy, reduction = "integrated.cca", split.by = "orig.ident")


#plot by reduction = "umap"
DimPlot(hthy, reduction = "umap", group.by = "orig.ident")
DimPlot(hthy, reduction = "umap", split.by = "orig.ident")

# find feature, save as pdf at 28in by 4in
FeaturePlot(hthy, reduction = "umap", features = "AIRE", split.by = "orig.ident")
FeaturePlot(hthy, reduction = "umap", features = "PRSS16", split.by = "orig.ident")
FeaturePlot(hthy, reduction = "umap", features = "PSMB11", split.by = "orig.ident")
FeaturePlot(hthy, reduction = "umap", features = "LY75", split.by = "orig.ident")
FeaturePlot(hthy, reduction = "umap", features = "COL4A6", split.by = "orig.ident")
FeaturePlot(hthy, reduction = "umap", features = "THY1", split.by = "orig.ident")
FeaturePlot(hthy, reduction = "umap", features = "ITGA6", split.by = "orig.ident")
FeaturePlot(hthy, reduction = "umap", features = "ITGB4", split.by = "orig.ident")
FeaturePlot(hthy, reduction = "umap", features = "BCAM", split.by = "orig.ident")




# test run 3D UMAP

# Interacive multimodal 3D UMAP plotting of scRNA sequencing datasets
# The following is a length of code generated to create nice 
# 3D UMAP plots of seurat v3.0.0-v3.1.1 objects utilizing the visualization 
# package plot_ly

# R v3.5.3 (x64 bit) and RStudio v1.2.1335 (x64 bit) were used for running this code :)

# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/

# Contributors (by their Github handles):
# @Dragonmasterx87 (Dept. of Cell Biology, UM)
# @msaadsadiq (Dept. of Electrical and Computer Engineering, UM)

# Install plot_ly
install.packages('plotly')

# Load plot_ly
library(plotly)

# Construct a dataframe using data from your pre-clustered Seurat v3.1.1 object
# Here 'seurat_clusters' is list of numeric cluster identities, you can find it here: hthy[["seurat_cluster"]], 
# or hthy$seurat_clusters, where 'hthy' is a Seurat object created with Seurat v3.1.1 (works for v3.0.0 as well)
hthy <- ThisIsWhateverhthyIsEvenIfItsIntegrated

# Re-run UMAPs that you have accurate calculations for all UMAP(s)
hthy <- RunUMAP(hthy, dims = 1:15, n.components = 3L)

# This is a manual method of extracting embeddings and is not needed
# as pointed out by user @sdinardo on 01142022 thank you! 
# Extract UMAP information from Seurat Object
# UMAP_1 <- hthy[["umap"]]@cell.embeddings[,1]
# UMAP_2 <- hthy[["umap"]]@cell.embeddings[,2]
# UMAP_3 <- hthy[["umap"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
head(Embeddings(object = hthy, reduction = "umap"))

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = hthy, vars = c("umap_1", "umap_2", "umap_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 17 clusters (0-16)
fig <- plot_ly(data = plot.data, 
               x = ~umap_1, y = ~umap_2, z = ~umap_3, 
               color = ~seurat_clusters, 
               colors = c("lightseagreen",
                          "gray50",
                          "darkgreen",
                          "red4",
                          "red",
                          "turquoise4",
                          "black",
                          "yellow4",
                          "royalblue1",
                          "lightcyan3",
                          "peachpuff3",
                          "khaki3",
                          "gray20",
                          "orange2",
                          "royalblue4",
                          "yellow3",
                          "gray80"),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 5, width=2), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


# Updates stemming from Issue #9 Having a fixed scale on axes while selecting particular clusters
# @rtoddler thanks for the suggestions!
# Before you plot, set the ranges of the axis you desire. This set axis range will be 
# present across all clusters, and plotly will not adjust for axis length anymore
# this axis length will persist even when selecting some clusters

# xaxis
axx <- list(
  nticks = 4,
  range = c(-10,10) #select range of xaxis
)

# yaxis
axy <- list(
  nticks = 4,
  range = c(-10,10) #select range of yaxis
)

#zaxis
axz <- list(
  nticks = 4,
  range = c(-10,10) #select range of zaxis
)

fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect
fig
fig_cube

# Say you wanto make a gene-expression 3D plot, where you can plot gene expression against a color scale
# Here using the same seurat object as above, we extract gene expression information for beta-actin 'ACTB'
# Here we concentrate on SCT normalized data, or log normalized RNA NOT raw counts.
# In addition if you want, you may look at normalised-RNA, SCT or integrated slots, to look at gene expression
# Setting your DefaultAssay() will inform R which assay to pick up expression data from.
DefaultAssay(object = hthy)
DefaultAssay(object = hthy) <- "RNA"
DefaultAssay(object = hthy) <- "integrated"
DefaultAssay(object = hthy) <- "SCT"

# create a dataframe
plot.data <- FetchData(object = hthy, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "ACTB"), slot = 'data')

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'changed' column of your dataframe
plot.data$changed <- ifelse(test = plot.data$ACTB <1, yes = plot.data$ACTB, no = 1)

# Add the label column, so that now the column has 'cellname-its expression value'
plot.data$label <- paste(rownames(plot.data)," - ", plot.data$ACTB, sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgreen', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5, width=2), 
        text=~label,
        hoverinfo="text"
)

# On running this code the HTML output should appear in RStudio. You can save the output as a
# HTML file. Once you have saved, just open the HTML file in any web browser (double click on the html- file
# and if asked select to open with any web browser like google chrome/safari/mozilla/explorer etc).
# It should be have all of the integrated features you saw in the RStudio output file.

########## #
########## #

# Alternative method as designed by @vertesy (Thanks for the suggestions!)
# create a dataframe
goi <- "TOP2A"
plotting.data <- FetchData(object = hthy, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression"=goi), slot = 'data')

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'Expression' column of your dataframe
# Cutoff <- 2
Cutoff <- quantile(plotting.data[,goi], probs = .95)
plotting.data$"ExprCutoff" <- ifelse(test = plotting.data[,goi] <Cutoff, yes = plotting.data[,goi], no = Cutoff)

# Add the label column, so that now the column has 'cellname-its expression value'
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data[,goi], sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plotting.data,
        # name = goi,
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~ExprCutoff, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgrey', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 1), 
        text=~label,
        hoverinfo="text"
) %>%layout(title=goi)

# Thank you for reading and using this code to further your scRNAseq analysis!
# If you liked it, dont forget to acknowledge, fork and star!
# Citation information is within the Readme, please dont forget to cite!
# Have a wonderful day!!





















############################################################################


# old script below for reference

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
saveRDS(hthy, file = "data/rds_objects/240205_parse_dataset_before_QC.rds")
hthy <- readRDS(file = "data/rds_objects/240205_parse_dataset_before_QC.rds")


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
dim(ht9)
dim(ht11)
dim(ht12)
dim(ht14)
dim(ht9_sub)
dim(ht11_sub)
dim(ht12_sub)
dim(ht14_sub)
# based on the dimensions above, determined that ht12 has 8032 cells, the lowest

# remove original objects
rm(ht9, ht11, ht12, ht14)


# load in sarah's reference data
# read in rds object with labels
reference <- readRDS(file = "data/rds_objects/meyer_snrna_seq_multiome_full_data_with_celltypes.rds")

# get the batch (sample) names from reference data set
batches <- unique(reference@meta.data$batch)

ref1 <- subset(x = reference, subset = batch == batches[1])
ref2 <- subset(x = reference, subset = batch == batches[2])
ref3 <- subset(x = reference, subset = batch == batches[3])

# merge all samples from parse and multiome dataset
hthy2 <- merge(x = ht9_sub, y = list(ht11_sub, 
                                     ht12_sub, ht14_sub, ref1, ref2, ref3, 
                                     ctec.cd49f, mtec.cd49f, ctec, mtec))

# merge only ht11 multiome and parse samples
ht11 <- merge(x = ht11_sub, y = ref2)

# standard analysis steps
hthy2 <- NormalizeData(hthy2)
hthy2 <- FindVariableFeatures(hthy2)
hthy2 <- ScaleData(hthy2)
hthy2 <- RunPCA(hthy2, features = VariableFeatures(object = hthy2))

hthy2 <- FindNeighbors(hthy2, dims = 1:30)
hthy2 <- FindClusters(hthy2, resolution = 0.5, 
                      cluster.name = "unintegrated_clusters")
hthy2 <- RunUMAP(hthy2, dims = 1:30, 
                 reduction.name = "umap.unintegrated")

# plot umap before integration
# group by parse samples ("sample")
DimPlot(hthy2, reduction = "umap.unintegrated", group.by = "sample")
# group by multiome samples ("batch")
DimPlot(hthy2, reduction = "umap.unintegrated", group.by = "batch")
# group by orig.ident
DimPlot(hthy2, reduction = "umap.unintegrated", group.by = "orig.ident")
DimPlot(hthy2, reduction = "umap.unintegrated", split.by = "orig.ident")


#save before integration
saveRDS(hthy2, file = "data/rds_objects/240908_before_integration.rds")
#hthy2 <- readRDS("data/rds_objects/240908_before_integration.rds")


# perform integration on the merged samples; wi = with integration
hthy2wi <- IntegrateLayers(object = hthy2, method = CCAIntegration, 
                           orig.reduction = "pca", 
                           new.reduction = "integrated.cca",
                           verbose = FALSE) # 240205 took 2h 38m 47s
                                            # 240906 took 1.5 hrs
                                            # 240908 took 2h 53s

# save object after integrating, before joining layers
saveRDS(hthy2wi, file = "data/rds_objects/240908_integrated_stem-tec.rds")
hthy2wi <- readRDS(file = "data/rds_objects/240908_integrated_stem-tec.rds")


# rejoin the layers after integration
hthy2wi[["RNA"]] <- JoinLayers(hthy2wi[["RNA"]])


# remove objects
rm(ht9_sub, ht11_sub, ht12_sub, ht14_sub)
rm(hthy, hthy2, ref1, ref2, ref3, reference)


# analyze post integration
hthy2wi <- FindNeighbors(hthy2wi, reduction = "integrated.cca", dims = 1:30)
hthy2wi <- FindClusters(hthy2wi, resolution = 0.5)
hthy2wi <- RunUMAP(hthy2wi, dims = 1:30, reduction = "integrated.cca")

# combine sample IDs
hthy2wi@meta.data <- hthy2wi@meta.data |>
  mutate(samples = ifelse(is.na(sample), batch, sample))

#plot graphs before combining sample IDs
DimPlot(hthy2wi, reduction = "umap", group.by = "sample")
DimPlot(hthy2wi, reduction = "umap", group.by = "batch")
DimPlot(hthy2wi, reduction = "umap", group.by = "samples")
DimPlot(hthy2wi, reduction = "umap", split.by = "sample")
DimPlot(hthy2wi, reduction = "umap", split.by = "batch")
DimPlot(hthy2wi, reduction = "umap", split.by = "samples")
DimPlot(hthy2wi, reduction = "umap", split.by = "orig.ident")

FeaturePlot(hthy2wi, reduction = "umap", features = "AIRE", split.by = "samples")
FeaturePlot(hthy2wi, reduction = "umap", features = "PRSS16", split.by = "samples")
FeaturePlot(hthy2wi, reduction = "umap", features = "PSMB11", split.by = "samples")
FeaturePlot(hthy2wi, reduction = "umap", features = "LY75", split.by = "samples")

FeaturePlot(hthy2wi, reduction = "umap", features = "AIRE", split.by = "sample")
FeaturePlot(hthy2wi, reduction = "umap", features = "PRSS16", split.by = "sample")
FeaturePlot(hthy2wi, reduction = "umap", features = "PSMB11", split.by = "sample")
FeaturePlot(hthy2wi, reduction = "umap", features = "LY75", split.by = "sample")
FeaturePlot(hthy2wi, reduction = "umap", features = "COL4A6", split.by = "sample")

FeaturePlot(hthy2wi, reduction = "umap", features = "AIRE", split.by = "batch")
FeaturePlot(hthy2wi, reduction = "umap", features = "PRSS16", split.by = "batch")
FeaturePlot(hthy2wi, reduction = "umap", features = "PSMB11", split.by = "batch")
FeaturePlot(hthy2wi, reduction = "umap", features = "LY75", split.by = "batch")
FeaturePlot(hthy2wi, reduction = "umap", features = "COL4A6", split.by = "batch")

FeaturePlot(hthy2wi, reduction = "umap", features = "AIRE", split.by = "orig.ident")
FeaturePlot(hthy2wi, reduction = "umap", features = "PRSS16", split.by = "orig.ident")
FeaturePlot(hthy2wi, reduction = "umap", features = "PSMB11", split.by = "orig.ident")
FeaturePlot(hthy2wi, reduction = "umap", features = "LY75", split.by = "orig.ident")
FeaturePlot(hthy2wi, reduction = "umap", features = "COL4A6", split.by = "orig.ident")
FeaturePlot(hthy2wi, reduction = "umap", features = "ITGA6", split.by = "orig.ident")
FeaturePlot(hthy2wi, reduction = "umap", features = "THY1", split.by = "orig.ident")
FeaturePlot(hthy2wi, reduction = "umap", features = "BCAM", split.by = "orig.ident")



# repeat analysis for merged ht11
ht11 <- NormalizeData(ht11)
ht11 <- FindVariableFeatures(ht11)
ht11 <- ScaleData(ht11)
ht11 <- RunPCA(ht11, features = VariableFeatures(object = ht11))

ht11 <- FindNeighbors(ht11, dims = 1:30)
ht11 <- FindClusters(ht11, resolution = 0.5, 
                      cluster.name = "unintegrated_clusters")
ht11 <- RunUMAP(ht11, dims = 1:30, 
                 reduction.name = "umap.unintegrated")

# plot umap before integration
# group by parse samples ("sample")
DimPlot(ht11, reduction = "umap.unintegrated", group.by = "sample")
# group by multiome samples ("batch")
DimPlot(ht11, reduction = "umap.unintegrated", group.by = "batch")
# group by orig.ident
DimPlot(ht11, reduction = "umap.unintegrated", group.by = "orig.ident")


# ht11 with integration
ht11wi <- IntegrateLayers(object = ht11, method = CCAIntegration, 
                           orig.reduction = "pca", 
                           new.reduction = "integrated.cca",
                           verbose = FALSE) # 240205 took 8m 18s

# save object after integrating and joining layers
saveRDS(ht11wi, file = "data/rds_objects/240205_ht11_integrated.rds")
ht11wi <- readRDS(file = "data/rds_objects/240205_ht11_integrated.rds")

# rejoin the layers after integration
ht11wi[["RNA"]] <- JoinLayers(ht11wi[["RNA"]])

# analyze post integration
ht11wi <- FindNeighbors(ht11wi, reduction = "integrated.cca", dims = 1:30)
ht11wi <- FindClusters(ht11wi, resolution = 0.5)
ht11wi <- RunUMAP(ht11wi, dims = 1:30, reduction = "integrated.cca")

# combine sample IDs
ht11wi@meta.data <- ht11wi@meta.data |>
  mutate(samples = ifelse(is.na(sample), batch, sample))

#plot graphs before combining sample IDs
DimPlot(ht11wi, reduction = "umap", group.by = "sample")
DimPlot(ht11wi, reduction = "umap", group.by = "batch")
DimPlot(ht11wi, reduction = "umap", group.by = "samples")
DimPlot(ht11wi, reduction = "umap", split.by = "samples")

FeaturePlot(ht11wi, reduction = "umap", features = "AIRE", split.by = "samples")
FeaturePlot(ht11wi, reduction = "umap", features = "PRSS16", split.by = "samples")
FeaturePlot(ht11wi, reduction = "umap", features = "PSMB11", split.by = "samples")
FeaturePlot(ht11wi, reduction = "umap", features = "LY75", split.by = "samples")


# save seurat object before trying label transfer
saveRDS(hthy2wi, file = "data/rds_objects/240205_all_7samples_integrated-pre_label_transfer.rds")
hthy2wi <- readRDS(file = "data/rds_objects/240205_all_7samples_integrated-pre_label_transfer.rds")
saveRDS(ht11wi, file = "data/rds_objects/240205_ht11_integrated-pre_label_transfer.rds")
ht11wi <- readRDS(file = "data/rds_objects/240205_ht11_integrated-pre_label_transfer.rds")



# label transfer for all 7 samples

tec.ref <- subset(hthy2wi, orig.ident == "SeuratProject")
tec.query <- subset(hthy2wi, orig.ident == "hthy")

tec.anchors <- FindTransferAnchors(reference = tec.ref, 
                                   query = tec.query,
                                   dims = 1:30, 
                                   reference.reduction = "pca")
predictions <- TransferData(anchorset = tec.anchors, 
                            refdata = tec.ref$celltype, 
                            dims = 1:30)
tec.query <- AddMetaData(tec.query, metadata = predictions)

# Visualization
DimPlot(tec.query, reduction = "umap")
DimPlot(tec.query, reduction = "umap", group.by = "samples")
DimPlot(tec.ref, reduction = "umap", group.by = "samples")
DimPlot(tec.query, reduction = "umap", split.by = "samples")
DimPlot(tec.ref, reduction = "umap", split.by = "samples")

# Unimodal UMAP Projection
tec.ref <- RunUMAP(tec.ref, dims = 1:30, 
                   reduction = "integrated.cca", 
                   return.model = TRUE)
tec.query <- MapQuery(anchorset = tec.anchors, 
                      reference = tec.ref, 
                      query = tec.query,
                      refdata = list(celltype = "celltype"), 
                      reference.reduction = "pca", reduction.model = "umap")

# plot predicted cell types
DimPlot(tec.query[,tec.query@meta.data$samples == "Ht9"], reduction = "umap", 
        group.by = "predicted.celltype", label = TRUE) + 
  NoLegend() +
  ggtitle("HT9")
DimPlot(tec.query[,tec.query@meta.data$samples == "HT11"], reduction = "umap", 
        group.by = "predicted.celltype", label = TRUE) +
  ggtitle("HT11") +
  NoLegend()
DimPlot(tec.query[,tec.query@meta.data$samples == "Ht12"], reduction = "umap", 
        group.by = "predicted.celltype", label = TRUE) +
  ggtitle("HT12") +
  NoLegend()
DimPlot(tec.query[,tec.query@meta.data$samples == "Ht14"], reduction = "umap", 
        group.by = "predicted.celltype", label = TRUE) +
  ggtitle("HT14") +
  NoLegend()
DimPlot(tec.ref, reduction = "umap", group.by = "celltype", label = TRUE) +
  ggtitle("Multiome Reference") +
  NoLegend()

# plot UMAP projection
p1 <- DimPlot(tec.ref, reduction = "umap", group.by = "celltype", 
              label = TRUE, label.size = 3,repel = TRUE,
              split.by = "samples") + 
  NoLegend() + 
  ggtitle("Reference annotations")
p2 <- DimPlot(tec.query, reduction = "ref.umap", 
              group.by = "predicted.celltype", 
              label = TRUE,
              label.size = 3, 
              repel = TRUE,
              split.by = "samples") + 
  NoLegend() + 
  ggtitle("Query transferred labels")
p1 + p2

p3 <- DimPlot(tec.ref, reduction = "umap", group.by = "celltype", 
              label = TRUE, label.size = 3,repel = TRUE) + 
  NoLegend() + 
  ggtitle("Reference annotations")
p4 <- DimPlot(tec.query, reduction = "ref.umap", 
              group.by = "predicted.celltype", 
              label = TRUE,
              label.size = 3, 
              repel = TRUE) + 
  NoLegend() + 
  ggtitle("Query transferred labels")
p3 + p4


tec.query.cellcounts <- table(tec.query$predicted.celltype, by=tec.query@meta.data$samples)
tec.ref.cellcounts <- table(tec.ref$celltype, by=tec.ref@meta.data$samples)

write.csv(tec.query.cellcounts, file = "240205_tec_query_tab.csv")
write.csv(tec.ref.cellcounts, file = "240205_tec_ref_tab.csv")




# label transfer for only ht11 sample

ht11.ref <- subset(ht11wi, orig.ident == "SeuratProject")
ht11.query <- subset(ht11wi, orig.ident == "hthy")

ht11.anchors <- FindTransferAnchors(reference = ht11.ref, 
                                   query = ht11.query,
                                   dims = 1:30, 
                                   reference.reduction = "pca")
predictions <- TransferData(anchorset = ht11.anchors, 
                            refdata = ht11.ref$celltype, 
                            dims = 1:30)
ht11.query <- AddMetaData(ht11.query, metadata = predictions)

# Visualization
DimPlot(ht11.query, reduction = "umap")
DimPlot(ht11..query, reduction = "umap", group.by = "samples")
DimPlot(ht11.ref, reduction = "umap", group.by = "samples")
DimPlot(ht11.query, reduction = "umap", split.by = "samples")
DimPlot(ht11.ref, reduction = "umap", split.by = "samples")

# Unimodal UMAP Projection
ht11.ref <- RunUMAP(ht11.ref, dims = 1:30, 
                   reduction = "integrated.cca", 
                   return.model = TRUE)
ht11.query <- MapQuery(anchorset = ht11.anchors, 
                      reference = ht11.ref, 
                      query = ht11.query,
                      refdata = list(celltype = "celltype"), 
                      reference.reduction = "pca", reduction.model = "umap")

# plot predicted cell types
DimPlot(ht11.query, reduction = "umap", 
        group.by = "predicted.celltype", split.by = "samples", label = TRUE) + 
  NoLegend()
DimPlot(ht11.query, reduction = "umap", 
        group.by = "predicted.celltype", label = TRUE) +
  ggtitle("HT11")

DimPlot(ht11.ref, reduction = "umap", group.by = "celltype", label = TRUE)

# plot UMAP projection
ht11.p1 <- DimPlot(ht11.ref, reduction = "umap", group.by = "celltype", 
              label = TRUE, label.size = 3,repel = TRUE,
              split.by = "samples") + 
  NoLegend() + 
  ggtitle("Reference annotations")
ht11.p2 <- DimPlot(ht11.query, reduction = "ref.umap", 
              group.by = "predicted.celltype", 
              label = TRUE,
              label.size = 3, 
              repel = TRUE,
              split.by = "samples") + 
  NoLegend() + 
  ggtitle("Query transferred labels")
ht11.p1 + ht11.p2

ht11.p3 <- DimPlot(ht11.ref, reduction = "umap", group.by = "celltype", 
              label = TRUE, label.size = 3,repel = TRUE) + 
  NoLegend() + 
  ggtitle("Reference annotations")
ht11.p4 <- DimPlot(ht11.query, reduction = "ref.umap", 
              group.by = "predicted.celltype", 
              label = TRUE,
              label.size = 3, 
              repel = TRUE) + 
  NoLegend() + 
  ggtitle("Query transferred labels")
ht11.p3 + ht11.p4


ht11.query.cellcounts <- table(ht11.query$predicted.celltype, by=ht11.query@meta.data$samples)
ht11.ref.cellcounts <- table(ht11.ref$celltype, by=ht11.ref@meta.data$samples)

write.csv(ht11.query.cellcounts, file = "240205_ht11_query_tab.csv")
write.csv(ht11.ref.cellcounts, file = "240205_ht11_ref_tab.csv")



# END OF SCRIPT 240205

# =================================

