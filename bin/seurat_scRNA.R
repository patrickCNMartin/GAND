#!/usr/bin/env Rscript

#-----------------------------------------------------------------------------#
# LIBRARIES
#-----------------------------------------------------------------------------#
library(future)
library(Matrix)
library(Seurat)
library(dplyr)
library(tidyr)
library(patchwork)
library(DESeq2)
library(ggplot2)
library(rmarkdown)
library(hdf5r)
library(ggpubr)
library(lme4)
library(emmeans)
library(argparser)
library(jsonlite)
set.seed(42)
#-----------------------------------------------------------------------------#
# ARGS & OPTIONS
#-----------------------------------------------------------------------------#
p <- arg_parser("Process input directory with manifest")

# Add arguments
p <- add_argument(p, "--input_dir", short = "-i", help = "Input directory", type = "character")
p <- add_argument(p, "--manifest", short = "-m", help = "Manifest file", type = "character")
p <- add_argument(p, "--ref_dir", short = "-r", help = "Input directory for referrence data", type = "character")
p <- add_argument(p, "--gene_sets", short = "-g",help = "Gene sets to check for enrichement", type = "character")
p <- add_argument(p, "--mut_genes", short = "-mg",help = "Mutally Excuslive Genes", type = "character")

# Parse arguments
argv <- parse_args(p)

print(argv)

# Assign to variables
input_dir <- argv$input_dir
manifest <- paste0(input_dir, "/", argv$manifest)
ref_dir <- argv$ref_dir
genes_sets <- fromJSON(argv$gene_sets)
mut_genes <- fromJSON(argv$mut_genes)

max_size <- 100000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# UTILS
#-----------------------------------------------------------------------------#
#' Load data into seurat objects
#' @param matrix_files vector{string} - location of mtx files
#' @param barcode_files vector{string} - location of barcode files
#' @param feature_files vector{string} - location of feature files
#' @param manifest data.frame - contains sample id and condition
#' @return list of seurat objects
load_data <- function(matrix_files,
                      barcode_files,
                      feature_files,
                      manifest) {
  seurat_objects <- vector("list", nrow(manifest))
  names(seurat_objects) <- manifest[,1]
  for (dat in seq_len(nrow(manifest))) {
    data_id <- manifest[dat, 1] # data set id
    mtx <- grep(data_id, matrix_files, value = TRUE)
    barcodes <- grep(data_id, barcode_files, value = TRUE)
    features <- grep(data_id, feature_files, value = TRUE)
    type <- rep(manifest[dat, 2],length(barcodes))
    counts <- ReadMtx(mtx = mtx,
                   features = features,
                   cells = barcodes)
    obj <- CreateSeuratObject(counts = counts)
    obj <- AddMetaData(obj,
                       metadata = type,
                       col.name = "condition")
    obj <- AddMetaData(obj,
                       metadata = data_id,
                       col.name = "sample")
    seurat_objects[[dat]] <- obj
  }
  return(seurat_objects)
}

#' Loading scRNA seq ref data with annotations
#' @param count_file path to counts
#' @param annotation_file path to cell annotatations.
#' @return Seurat object with pre-processed data 
load_ref <- function(count_file,
                     annotation_file,
                     ref_tag,
                     dims = 1:30) {
  # Loading counts
  ref_counts <- read.table(count_file,
                          header = TRUE,
                          row.names = 1,
                          sep = "\t",
                          check.names = FALSE)
  ref_counts <- CreateSeuratObject(counts = ref_counts,
                                   project =  ref_tag)
  ref_counts <- ref_counts %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>% 
    RunUMAP(dims = dims) # why 22 dims?
  
  # Loading annotations
  annot <- read.delim(annotation_file,
                     header = TRUE,
                     check.names = FALSE)
  rownames(annot) <- annot$CellID
  annot$CellID <- NULL
  rownames(annot) <- gsub("e14\\.", "e14-", rownames(annot))
  rownames(annot) <- gsub("e14-WT9\\.", "e14-WT9-", rownames(annot))
  rownames(annot) <- gsub("e14-WT8\\.", "e14-WT8-", rownames(annot))
  annot <- annot[match(colnames(ref_counts), rownames(annot)),]
  # Add annotations
  ref_counts <- AddMetaData(ref_counts, metadata = annot, col.name = "cell_type")
  return(ref_counts)
}

#' QC seurat data
#' @param seurat_object seurat object containing counts
#' @param pattern regex string - define mitochondrial gene patterns
#' @param feature_range vector{int} - min and max number of genes to retain per cell
#' @param percent.mt int - max percentage of mitochondrial genes allowed per cell
#' @param nfeatures int - number of variable features 
#' @param npcs int - number of prinical componenents 
#' @return seurat object
seurat_qc <- function(seurat_object,
                      pattern = "^mt-",
                      feature_range = c(0, 10000),
                      percent_mt = 10,
                      nfeatures = 2000,
                      npcs = 30,
                      skip_subset = FALSE) {
  if (!skip_subset) {
    # Setting Mito percentage
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
                                                        pattern = pattern)
    # Subsetting object
    min_features <- min(feature_range)
    max_features <- max(feature_range)
    seurat_object <- subset(seurat_object,
                            subset = nFeature_RNA > min_features &
                            nFeature_RNA < max_features &
                            percent.mt < percent_mt)
  }
  
  # Process data
  seurat_object <- seurat_object %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = nfeatures) %>%
    ScaleData() %>%
    RunPCA()
  # Return Object
  return(seurat_object)
}
#' perform clustering and UMAP projections on seurat objecy
#' @param seurat_object pre-processd seurat object
#' @param dim_n int - number of PCA dims to use for clustering
#' @param reduction string - which reduced dim space to for clustering
#' @param resolution numeric - louvain clustering resolution
#' @param cluster_name string - name tag for clusters
#' @return seurat object
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

#' integrate seurat object list
#' @param seurat_list list - list of QC'ed seurat objects
#' @param method string - integration method
#' @param reduction string - dim reduction method to use for integration
#' @param integration_tag string - tag used for new inetgrated dim reduction assay
#' @param dim_n int - number of PCs to use during integration and clustering
#' @param resolution numeric - louvain clustering resolution
#' @param cluster_name string - tag used for integrated clusters
#' @return integrated seurat object
integrate_list <- function(seurat_list,
                           method = "RPCAIntegration",
                           reduction = "pca",
                           integration_tag = "integrated",
                           dim_n = 30,
                           resolution = 0.4,
                           cluster_name = "integrated_cluster_"){
  seurat_merged <- merge(seurat_list[[1]], seurat_list[-1])
  seurat_merged <- seurat_qc(seurat_merged, skip_subset = TRUE)
  seurat_merged <- IntegrateLayers(object = seurat_merged,
                                   method = method,
                                   orig.reduction = reduction,
                                   new.reduction = integration_tag,
                                   verbose = FALSE)
  #seurat_merged <- JoinLayers(seurat_merged)
  seurat_merged <- seurat_clusters(seurat_merged,
                                  dim_n = dim_n,
                                  reduction = integration_tag,
                                  resolution = resolution,
                                  cluster_name = cluster_name)
  # Combinding labels 
  combinded_type <- seurat_merged@meta.data$condition
  combinded_type[grep("WT",combinded_type)] <- "WT"
  combinded_type[grep("KO",combinded_type)] <- "KO"
  combinded_type[grep("Het",combinded_type)] <- "Het"
  seurat_merged <- AddMetaData(seurat_merged,
                     metadata = combinded_type,
                     col.name = "condition_aggregated")
  return(seurat_merged)
}

transfer_labels <- function(scrna,
                            ref,
                            dims = 1:30,
                            reduction = "pca") {
  anchors <- FindTransferAnchors(reference = ref,
                                 query = scrna,
                                 dims = dims,
                                 reference.reduction = reduction)
  #TransferData 
  predictions <- TransferData(anchorset = anchors,
                              refdata = ref$cell_type,
                              dims = dims)
  scrna <- AddMetaData(scrna,
                       metadata = predictions)
  return(scrna)

}

#-----------------------------------------------------------------------------#
# DATA LOADING
#-----------------------------------------------------------------------------#
## scRNA data
manifest <- read.csv(manifest, header = FALSE, sep = " ")
mtx_files <- list.files(path = input_dir,
                        pattern = "matrix.mtx.gz",
                        full.names = TRUE)
feat_files <- list.files(path = input_dir,
                         pattern = "features.tsv.gz",
                         full.names = TRUE)
barcode_files <- list.files(path = input_dir,
                            pattern = "barcodes.tsv.gz",
                            full.names = TRUE)


seurat_list <- load_data(matrix_files = mtx_files,
                        barcode_files = barcode_files,
                        feature_files = feat_files,
                        manifest = manifest)

## Reference data
count_file <- list.files(path = ref_dir,
                        pattern = "_combined_matrix.txt.gz",
                        full.names = TRUE)
ref_tag  <- list.files(path = ref_dir,
                        pattern = "_combined_matrix.txt.gz",
                        full.names = FALSE)
ref_tag <- gsub("_combined_matrix.txt.gz", "", ref_tag)
ref_annot <- list.files(path = ref_dir,
                        pattern = "_combined_matrix_ClusterAnnotations.txt.gz",
                        full.names = TRUE)

ref_data <- load_ref(count_file,
                    ref_annot,
                    ref_tag)

#-----------------------------------------------------------------------------#
# QC and save
#-----------------------------------------------------------------------------#
seurat_list <- lapply(seurat_list,
                      FUN = seurat_qc,
                      pattern = "^mt-",
                      feature_range = c(0, 10000),
                      percent_mt = 10,
                      nfeatures = 2000,
                      npcs = 30)

seurat_list <- lapply(seurat_list,
                      FUN = seurat_clusters,
                      dim_n = 30,
                      reduction = "pca",
                      resolution = 0.4,
                      cluster_name = NULL)
saveRDS(seurat_list, file = "GAND_preprocessed.rds")

#-----------------------------------------------------------------------------#
# Integration
#-----------------------------------------------------------------------------#
seurat_integrated <- integrate_list(seurat_list,
                                    method = "RPCAIntegration",
                                    reduction = "pca",
                                    integration_tag = "integrated",
                                    dim_n = 30,
                                    resolution = 0.4,
                                    cluster_name = "integrated_cluster_")

saveRDS(seurat_integrated, file = "GAND_seurat_integrated.rds")
#-----------------------------------------------------------------------------#
# Transfer labels
#-----------------------------------------------------------------------------#
seurat_annotated <- transfer_labels(seurat_integrated,
                                   ref_data)
saveRDS(seurat_annotated, file = "GAND_seurat_annotated.rds")
#-----------------------------------------------------------------------------#
# DONE
#-----------------------------------------------------------------------------#

