##Before we can analyse the data, some pre processing is necessary
##the read of bad quality have been removed, but for examples, some cells have too low "features" (i.e. genes)

library(Seurat)

pbmc.data = Read10X(data.dir = "../input_data/10X_CellRanger/filtered_feature_bc_matrix/")
#Create the Seurat Object
#TODO explain what the difference between pbmc.data and pbmc
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pouloup", min.cells = 3, min.features = 20)
#In creating the Seurat object, a first (careful) filter is already made :
#We only keep cells that have at least 20 features (i.e. different genes)
#Similarly, we only keep the

print("poutoup")

##Check and compute important features
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head (pbmc@meta.data, 5)

##TODO Show off some other options possible


png("features_plot.png")
VlnPlot(pbmc, c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()


#Okay so it's ... bimodal ? let's check wether it's correlated (it should be)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

png("scatter_plots.png", height = 20, width = 20, units="cm", res=100)
plot1 + plot2
dev.off()

#The bimodal stuff does not look like it is from two cells in the same bubble (cause they would share features, so total number of features should not be twice as high)

pbmc <- subset(pbmc, subset  = nFeature_RNA > 200 &  nFeature_RNA < 6000 & percent.mt <16)
png("features_plot_after_filter.png")
VlnPlot(pbmc, c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()
