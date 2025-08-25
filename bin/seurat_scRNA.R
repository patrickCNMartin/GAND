#-----------------------------------------------------------------------------#
# LIBRARIES
#-----------------------------------------------------------------------------#
library(remotes)
library(future)
library(Matrix)
library(Seurat)
library(dplyr)
library(patchwork)
library(DESeq2)
library(presto)
library(hdf5r)
library(ggpubr)
set.seed(42)

#-----------------------------------------------------------------------------#
# UTILS
#-----------------------------------------------------------------------------#

seurat_qc <- function(seurat_object,
                      pattern = "^mt-",
                      feature_range = c(0, 10000),
                      percent.mt = 10,
                      nfeatures = 2000) {
  # Setting Mito percentage
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
                                                        pattern = pattern)
  # Subsetting object
  min_features <- min(feature_range)
  max_features <- max(feature_range)
  seurat_object <- subset(seurat_object,
                          subset = nFeature_RNA > min_features &
                            nFeature_RNA < max_features &
                            percent.mt < percent.mt)
  # Process data
  seurat_object <- seurat_object %>%
    NormalizeData(seurat_object) %>%
    FindVariableFeatures(nfeatures = nfeatures) %>%
    ScaleData() %>%
    RunPCA()
  # Return Object
  return(seurat_object)
}


seurat_clusters <- function(seurat_object,
                            dim_n = 30,
                            reduction = "pca",
                            resolution = 0.4,
                            cluster_name = NULL) {
  if (!is.null(cluster_name)) {
    cluster_name <- paste0(cluster_name, resolution)
  }
  seurat_object <- seurat_object %>%
    FindNeighbors(dims = 1:dim_n, reduction = reduction) %>%
    FindClusters(resolution = resolution, cluster.name = cluster_name) %>%
    RunUMAP(dims = 1:dim_n)
  return(seurat_object)
}

#-----------------------------------------------------------------------------#
# ARGS
#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
input_dir <- args[1]
manifest <- args[2]
#-----------------------------------------------------------------------------#
# DATA LOADING
#-----------------------------------------------------------------------------#

mtx_files <- list.files(path = input_dir, pattern = "matrix.mtx.gz", full.names = TRUE)
feat_files <- list.files(path = input_dir, pattern = "features.tsv.gz", full.names = TRUE)
cells_files <- list.files(path = input_dir, pattern = "barcodes.tsv.gz", full.names = TRUE)
sample_names <- gsub("_matrix.mtx.gz", "", mtx_files)


seurat_list <- vector(mode = "list", length = length(mtx_files))
names(seurat_list) <- paste0(sample_names, "_seurat_object")


for (i in seq_along(mtx_files)) {
  expression_matrix <- ReadMtx(mtx = mtx_files[i],
                               features = feat_files[i],
                               cells = cells_files[i])
  seurat_object <- CreateSeuratObject(counts = expression_matrix)
  
  seurat_object <- seurat_qc(seurat_object)
  
  sample <- rep(sample_names[i], length(cells_files[i]))
  seurat_object <- AddMetaData(seurat_object,
                               metadata = sample,
                               col.name = "sample")
  seurat_list[[i]] <- seurat_object
}






# We will use the for loop here - below I show some otherways you
# can do loops in a cleaner way if you don't need to many things

plot_list <- vector("list", length(seurat_list))
for(so in seq_along(seurat_list)) {
  # As it is the case with assign - avoid using get() - powerful
  # but not required here. 
  # SO <- get(seurat_objects[j])
  # In this case, you could have use SO <- seurat_object[j]
  # We are using a list so we can simple use SO <- seurat_list[[j]]
  # We already ran PCA in the seurat_qc function so we 
  # Note that here the Elbowplot returns a ggplot object
  # Default number of PCs in Seurat v5 is 50 
  plot_list[[so]] <- ElbowPlot(seurat_list[[so]], ndims = 50)
}

# Using ggpurb to plot all pannels
shape <- ceiling(sqrt(length(plot_list)))
plot_list <- ggarrange(plotlist = plot_list,
                       ncol = shape,
                       nrow = shape)
print(plot_list)


# Step 8 ---- Step 8. Cluster the cells


# This section it could be easier to save plots
# You don't need to rep the PCA dim - you don't need to loop
# over that variable. 
# For is it could be 
for (k in seq_along(seurat_list)){
  p <- PCHeatmap(SO, dims = 1:20, cells = 500, balanced = TRUE, ncol = 4)
  print(p)
}




# I have put the clustering sepperately since it could change the number
# of dims you will choose
# As a side note - consider using apply loops 
# here is an example 

seurat_list <- lapply(seurat_list, seurat_cluster, dim_n = 20, resolution = 0.4)


umap_list <- vector("list", length(seurat_list))
for(so in seq_along(seurat_list)) {
  umap_list[[so]] <- DimPlot(seurat_list[[so]], reduction = "umap")
}

# Using ggpurb to plot all pannels
shape <- ceiling(sqrt(length(umap_list)))
umap_list <- ggarrange(plotlist = umap_list,
                       ncol = shape,
                       nrow = shape)
print(umap_list)

# I would recommend using rds files instead of Rdata/Rda 
# Ironically this is where the get function would come in handy if you were
# to reload the object. Rdata works a bit strangely when reloading objects
# save RDS will save a more effcient version of the data and when you reload it
# using readRDS("file.rds"), you can assign it to a new variable. 
# this store the entire list.
saveRDS(seurat_list, file = paste0(file_dir,"Seurat_Objects_for_Annotation.rds"))






# This is where you using get is useful and required
# However - we don't need to worry with all of this
# Let's imagine that you loading your data from the rds
# Now everything is already in a list and you don't need to
# search through the environment 
seurat_list <- readRDS(paste0(file_dir,"Seurat_Objects_for_Annotation.rds"))

## First we merge all objects 
seurat_merged <- merge(seurat_list)

## We then proceed with a joint analysis - this is running the same as on
## individual data setsseurat_list
seurat_merged <- seurat_qc(seurat_merged)

## Now we integrate the data 
## NOTE: I am using RPCAIntegration but seurat offers
## A wide variety of other methods
seurat_merged <- IntegrateLayers(object = seurat_merged,
                                 method = "RPCAIntegration",
                                 orig.reduction = "pca",
                                 new.reduction = "integrated.rpca",
                                 verbose = FALSE)

seurat_merged <- seurat_cluster(seurat_merged,
                            dim_n = 30,
                            reduction = "integrated.rpca",
                            resolution = 0.4,
                            cluster_name = "integrated_cluster_")

#seurat_merged[["RNA"]] <- JoinLayers(seurat_merged)

p1 <- DimPlot(seurat_merged, reduction = "umap", group.by = "sample")
p2 <- DimPlot(seurat_merged, reduction = "umap", group.by = "integrated_cluster_0.4")

print(p1 + p2)
