##Before we can analyse the data, some pre processing is necessary
##the read of bad quality have been removed, but for examples, some cells have too low "features" (i.e. genes)
#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
library(Seurat)
require(data.table)
require(ggplot2)
pbmc.data = Read10X(data.dir = "../input_data/seurat_tutorial/filtered_gene_bc_matrices/hg19")
#Create the Seurat Object
#TODO explain what the difference between pbmc.data and pbmc
#Save and reload pbmc object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pouloup", min.cells = 3, min.features = 20)
#In creating the Seurat object, a first (careful) filter is already made :
#We only keep cells that have at least 20 features (i.e. different genes)
#Similarly, we only keep the

print("poutoup")

##TODO A WHOLE THING on doublet detections and principle

##Check and compute important features
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head (pbmc@meta.data, 5)

##TODO Show off some other options possible


png("../img/features_plot.png")
VlnPlot(pbmc, c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()


#Okay so it's ... bimodal ? let's check wether it's correlated (it should be)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

png("../img/scatter_plots.png", height = 20, width = 20, units="cm", res=100)
plot1 + plot2
dev.off()

#The bimodal stuff does not look like it is from two cells in the same bubble (cause they would share features, so total number of features should not be twice as high)
#But it does very strange stuff with the PCs, so we'll cut it out

pbmc <- subset(pbmc, subset  = nFeature_RNA > 200 &  nFeature_RNA < 2500 & percent.mt <5)
png("../img/features_plot_after_filter.png")
VlnPlot(pbmc, c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()



##Data normalisation
##TODO Speak of various normalisation methods and their various assumptions and practical uses.
##Their little talk is already pretty good, but a page on "why do we normalise data" could be useful

pbmc <- NormalizeData(pbmc)



##Highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png("../img/highly_variable_features.png", width = 30, height = 20, units = "cm", res = 100)
plot1 + plot2
dev.off()


##Scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#TODO write a paragraph about what is scaling the data
#TODO Write a paragraph for the interests and problems with scaling the data
##TODO : add options to regress variations from unwanted sources (mt.DNA, Batch, GC ?)


##PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#TODO : Explanation of a PCA
#PCA only run on variable features that have been scaled, so here, only the first 2000 !
#Why only variables features
#Ways to look at the PCA results
#TODO : exxplain that Seurat uses PCA results for a bunch of other things, so it has to be computed
print(pbmc[["pca"]], dims = 1:5, nfeatures = 4)

png("../img/genes_participations_to_PCs.png", width = 30, height = 20, units = "cm", res = 100)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
dev.off()

png("../img/pca.png", width = 20, height = 20, units = "cm", res = 100)
DimPlot(pbmc, reduction = "pca")
dev.off()

png("../img/PCA_by_cells.png", width = 30, height = 20, units="cm", res=100)
DimHeatmap(pbmc, dims = 1, cells = 100, balanced = TRUE)
dev.off()


##Determine the dimensionality of the dataset (i.e. clusterisation)
##TODO Why it's needed
##Descirption of the method
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

png("../img/JackStrawPlot.png", width = 20, height = 20, units = "cm", res = 100)
JackStrawPlot(pbmc, dims = 1:20)
dev.off()


png("../img/Elbowplot.png", width = 20, height = 20, units = "cm", res = 100)
ElbowPlot(pbmc)
dev.off()


##Clustering
#distance is computed based on PCs
#TODO : curse of dimensionality with the dimensionality reduction thing ?
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)


##UMAP

pbmc <- RunUMAP(pbmc, dims = 1:10) #Dims here are still PCA decided
png("../img/UMAP.png", width = 20, height = 20, units = "cm", res = 100)
DimPlot(pbmc, reduction = "umap")
dev.off()
