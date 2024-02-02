# Loading libraries and setting location paths
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)

# remove all objects (variables, functions, etc.) from the current working 
# environment or session.
rm(list = ls())


# install SeuratData package
# devtools::install_github('satijalab/seurat-data')

# load SeuratData package
library(SeuratData)

# install data
# InstallData("panc8")
# InstallData("ifnb")

# load data
ifnb <- LoadData("ifnb")

# split the RNA measurements into two layers one for control cells, one for
# stimulated cells
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f=ifnb$stim)


# perform analysis without integration

# run standard anlaysis workflow
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)

ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))



# perform integration

ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")

# Visualization
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))
# To visualize the two conditions side-by-side, we can use the split.by argument
# to show each condition colored by cluster.
DimPlot(ifnb, reduction = "umap", split.by = "stim")


# To identify canonical cell type marker genes that are conserved across conditions, 
# we provide the FindConservedMarkers() function. This function performs differential 
# gene expression testing for each dataset/group and combines the p-values using 
# meta-analysis methods from the MetaDE R package. For example, we can calculated 
# the genes that are conserved markers irrespective of stimulation condition in 
# cluster 6 (NK cells).

# identify canonical cell type marker genes that are conserved across conditions
Idents(ifnb) <- "seurat_annotations"
nk.markers <- FindConservedMarkers(ifnb, ident.1 = "NK", grouping.var = "stim", verbose = FALSE)
head(nk.markers)

# install a faster package for wilcoxon rank sum test
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')


# You can perform these same analysis on the unsupervised clustering results 
# (stored in seurat_clusters), and use these conserved markers to annotate cell 
# types in your dataset.

# The DotPlot() function with the split.by parameter can be useful for viewing 
# conserved cell type markers across conditions, showing both the expression level
# and the percentage of cells in a cluster expressing any given gene. Here we plot
# 2-3 strong marker genes for each of our 14 clusters.

# NEEDS TO BE FIXED AND SET ORDER CORRECTLY
# NEEDS TO BE FIXED AND SET ORDER CORRECTLY
Idents(ifnb) <- factor(Idents(ifnb), levels = c("pDC", "Eryth", "Mk", "DC", 
                                                "CD14 Mono", "CD16 Mono",
                                                "B Activated", "B", "CD8 T", 
                                                "NK", "T activated", 
                                                "CD4 Naive T", "CD4 Memory T"))

markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY",
                     "NKG7", "CCL5","CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1",
                     "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1","GPR183", 
                     "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", 
                     "PRSS57")

DotPlot(ifnb, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, 
        split.by = "stim") + RotatedAxis()

FeaturePlot(ifnb, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", 
            max.cutoff = 3, cols = c("grey", "red"), reduction = "umap")


# Identify differential expressed genes across conditions

# Now that we’ve aligned the stimulated and control cells, we can start to do 
# comparative analyses and look at the differences induced by stimulation.

# We can aggregate cells of a similar type and condition together to create 
# “pseudobulk” profiles using the AggregateExpression command. As an initial 
# exploratory analysis, we can compare pseudobulk profiles of two cell types 
# (naive CD4 T cells, and CD14 monocytes), and compare their gene expression 
# profiles before and after stimulation. We highlight genes that exhibit dramatic 
# responses to interferon stimulation. As you can see, many of the same genes are 
# upregulated in both of these cell types and likely represent a conserved 
# interferon response pathway, though CD14 monocytes exhibit a stronger 
# transcriptional response.

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

aggregate_ifnb <- AggregateExpression(ifnb, group.by = c("seurat_annotations", 
                                                "stim"), return.seurat = TRUE)
genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", 
                   "CXCL10", "CCL8")

p1 <- CellScatter(aggregate_ifnb, "CD14 Mono_CTRL", "CD14 Mono_STIM", 
                  highlight = genes.to.label)
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

p3 <- CellScatter(aggregate_ifnb, "CD4 Naive T_CTRL", "CD4 Naive T_STIM", 
                  highlight = genes.to.label)
p4 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)

p2 + p4


# We can now ask what genes change in different conditions for cells of the same 
# type. First, we create a column in the meta.data slot to hold both the cell type 
# and stimulation information and switch the current ident to that column. Then we 
# use FindMarkers() to find the genes that are different between stimulated and 
# control B cells. Notice that many of the top genes that show up here are the 
# same as the ones we plotted earlier as core interferon response genes. 
# Additionally, genes like CXCL10 which we saw were specific to monocyte and B 
# cell interferon response show up as highly significant in this list as well.

# Please note that p-values obtained from this analysis should be interpreted with
# caution, as these tests treat each cell as an independent replicate, and ignore 
# inherent correlations between cells originating from the same sample. As 
# discussed here, DE tests across multiple conditions should expressly utilize 
# multiple samples/replicates, and can be performed after aggregating 
# (‘pseudobulking’) cells from the same sample and subpopulation together. We do
# not perform this analysis here, as there is a single replicate in the data, but 
# please see our vignette comparing healthy and diabetic samples as an example for 
# how to perform DE analysis across conditions.

ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
b.interferon.response <- FindMarkers(ifnb, ident.1 = "B_STIM", ident.2 = 
                                       "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)


# Another useful way to visualize these changes in gene expression is with the 
# split.by option to the FeaturePlot() or VlnPlot() function. This will display 
# FeaturePlots of the list of given genes, split by a grouping variable 
# (stimulation condition here). Genes such as CD3D and GNLY are canonical cell 
# type markers (for T cells and NK/CD8 T cells) that are virtually unaffected by 
# interferon stimulation and display similar gene expression patterns in the 
# control and stimulated group. IFI6 and ISG15, on the other hand, are core 
# interferon response genes and are upregulated accordingly in all cell types. 
# Finally, CD14 and CXCL10 are genes that show a cell type specific interferon 
# response. CD14 expression decreases after stimulation in CD14 monocytes, which 
# could lead to misclassification in a supervised analysis framework, underscoring
# the value of integrated analysis. CXCL10 shows a distinct upregulation in 
# monocytes and B cells after interferon stimulation but not in other cell types.

FeaturePlot(ifnb, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", 
            max.cutoff = 3, cols = c("grey","red"), reduction = "umap")

plots <- VlnPlot(ifnb, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim",
                 group.by = "seurat_annotations",pt.size = 0, combine = FALSE)

wrap_plots(plots = plots, ncol = 1)


# Perform integration with SCTransform-normalized datasets

# As an alternative to log-normalization, Seurat also includes support for 
# preprocessing of scRNA-seq using the sctransform workflow. The IntegrateLayers
# function also supports SCTransform-normalized data, by setting the 
# normalization.method parameter, as shown below.

ifnb <- LoadData("ifnb")

# split datasets and process without integration
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb <- SCTransform(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- RunUMAP(ifnb, dims = 1:30)
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))

# integrate datasets
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, 
                        normalization.method = "SCT", verbose = F)
ifnb <- FindNeighbors(ifnb, reduction = "integrated.dr", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 0.6)

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.dr")
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))

# perform differential expression
ifnb <- PrepSCTFindMarkers(ifnb)
ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
b.interferon.response <- FindMarkers(ifnb, ident.1 = "B_STIM",
                                     ident.2 = "B_CTRL", verbose = FALSE)