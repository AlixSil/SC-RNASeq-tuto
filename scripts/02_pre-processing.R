##Before we can analyse the data, some pre processing is necessary
##the read of bad quality have been removed, but for examples, some cells have too low "features" (i.e. genes)

library(Seurat)

NSCLC_DTC.data = Read10X(data.dir = "../input_data/10X_CellRanger/sample_feature_bc_matrix/")
str(NSCLC_DTC.data)
print("poutoup")
