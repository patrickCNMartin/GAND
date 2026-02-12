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
library(rhdf5)
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
p <- add_argument(p, "--manifest", short = "-m", help = "Path to manifest file", type = "character")
p <- add_argument(p, "--ref_dir", short = "-r", help = "Input directory for referrence data", type = "character")
p <- add_argument(p, "--number_pcs", short = "-n", help = "Number of Principal Components", type = "numeric")
p <- add_argument(p, "--min_features", short = "-minf", help = "Minimum number of features per cell", type = "numeric")
p <- add_argument(p, "--max_features", short = "-maxf", help = "Maximum number of features per cell", type = "numeric")
p <- add_argument(p, "--percent_mt", short = "-pmt", help = "Percentage of Mitochondrial RNA allowed", type = "numeric")
p <- add_argument(p, "--n_var_features", short = "-nf", help = "Number of variable Features", type = "numeric")
p <- add_argument(p, "--cluster_resolution", short = "-nf", help = "Louvain Clustering Resolution", type = "numeric")
p <- add_argument(p, "--integration_tag", short = "-it", help = "Tag used for integration", type = "character")
p <- add_argument(p, "--integration_method", short = "-im", help = "Integration Method", type = "character")

# Parse arguments
argv <- parse_args(p)

# Parse arguments
args <- parse_args(p)

# Convert to variables
input_dir <- args$input_dir
manifest <- args$manifest
ref_dir <- args$ref_dir
number_pcs <- args$number_pcs
min_features <- args$min_features
max_features <- args$max_features
percent_mt <- args$percent_mt
n_var_features <- args$n_var_features
cluster_resolution <- args$cluster_resolution
integration_tag <- args$integration_tag
integration_method <- args$integration_method


max_size <- 100000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# UTILS
#-----------------------------------------------------------------------------#
#' Load Data into Seurat Objects
#' 
#' Reads MTX files based on a manifest and creates a list of Seurat objects 
#' with sample and condition metadata.
#' 
#' @param matrix_files Character vector of paths to .mtx files.
#' @param barcode_files Character vector of paths to barcode files.
#' @param feature_files Character vector of paths to feature/gene files.
#' @param manifest Data frame containing sample IDs in column 1 and conditions in column 2.
#'
#' @return A named list of Seurat objects.
#' @export
load_data <- function(matrix_files,
                      barcode_files,
                      feature_files,
                      manifest) {
  seurat_objects <- vector("list", nrow(manifest))
  names(seurat_objects) <- manifest[,1]
  for (dat in seq_len(nrow(manifest))) {
    data_id <- manifest[dat, 1] 
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

#' Load Reference scRNA-seq Data with Annotations
#' 
#' Loads a count matrix and annotation file to create a processed reference 
#' Seurat object, including normalization, PCA, and UMAP.
#' 
#' @param count_file Path to the tab-separated count matrix file.
#' @param annotation_file Path to the tab-separated cell annotation file.
#' @param ref_tag Character string for the project name.
#' @param dims Integer vector specifying dimensions to use for UMAP.
#'
#' @return A processed Seurat object with a 'cell_type' metadata column.
#' @export
load_ref <- function(count_file,
                     annotation_file,
                     ref_tag,
                     dims = 1:30) {
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
    RunUMAP(dims = dims)
  
  annot <- read.delim(annotation_file,
                     header = TRUE,
                     check.names = FALSE)
  rownames(annot) <- annot$CellID
  annot$CellID <- NULL
  
  rownames(annot) <- gsub("e14\\.", "e14-", rownames(annot))
  rownames(annot) <- gsub("e14-WT9\\.", "e14-WT9-", rownames(annot))
  rownames(annot) <- gsub("e14-WT8\\.", "e14-WT8-", rownames(annot))
  annot <- annot[match(colnames(ref_counts), rownames(annot)),]
  
  ref_counts <- AddMetaData(ref_counts, metadata = annot, col.name = "cell_type")
  return(ref_counts)
}

#' Quality Control and Processing of Seurat Objects
#' 
#' Filters cells based on mitochondrial percentage and feature counts, 
#' then performs standard normalization and PCA.
#' 
#' @param seurat_object A Seurat object.
#' @param pattern Regex pattern to identify mitochondrial genes (default: "^mt-").
#' @param min_feature Minimum number of features required to keep a cell.
#' @param max_feature Maximum number of features allowed per cell.
#' @param percent_mt Maximum allowed mitochondrial percentage per cell.
#' @param nfeatures Number of variable features to find.
#' @param npcs Number of principal components to compute.
#'
#' @return A filtered and PCA-transformed Seurat object.
#' @export
seurat_qc <- function(seurat_object,
                      pattern = "^mt-",
                      min_feature = 0,
                      max_feature = 10000,
                      percent_mt = 10,
                      nfeatures = 2000,
                      npcs = 30) {
  
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
                                                        pattern = pattern)
    
    seurat_object <- subset(seurat_object,
                            subset = nFeature_RNA > min_feature &
                            nFeature_RNA < max_feature &
                            percent.mt < percent_mt)
  
  seurat_object <- seurat_object %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = nfeatures) %>%
    ScaleData() %>%
    RunPCA(npcs = npcs)
    
  return(seurat_object)
}

#' Perform Clustering and UMAP Projection
#' 
#' Finds neighbors, calculates clusters using Louvain algorithm, and generates UMAP.
#' 
#' @param seurat_object A pre-processed Seurat object.
#' @param dim_n Number of dimensions to use for neighbors and UMAP.
#' @param reduction Reduction space to use (e.g., "pca" or "integrated").
#' @param resolution Value for Louvain clustering granularity.
#' @param cluster_name Optional string tag for the metadata column name.
#'
#' @return A Seurat object with cluster assignments and UMAP coordinates.
#' @export
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

#' Integrate a List of Seurat Objects
#' 
#' Merges a list of objects, performs integration (default RPCA), and 
#' creates aggregated condition labels (WT, Het, KO).
#' 
#' @param seurat_list A list of Seurat objects.
#' @param method Integration method name (e.g., "RPCAIntegration").
#' @param reduction The original reduction to use for integration.
#' @param integration_tag Name for the new integrated reduction.
#' @param dim_n Number of PCs for integration and clustering.
#' @param resolution Louvain clustering resolution.
#' @param cluster_name Tag prefix for the integrated clusters.
#'
#' @return An integrated and clustered Seurat object.
#' @export
integrate_list <- function(seurat_list,
                           method = "RPCAIntegration",
                           reduction = "pca",
                           integration_tag = "integrated",
                           dim_n = 30,
                           resolution = 0.4,
                           cluster_name = "integrated_cluster_"){
  seurat_merged <- merge(seurat_list[[1]], seurat_list[-1])
  
  # Note: ensure seurat_qc handles a merged object correctly
  seurat_merged <- seurat_qc(seurat_merged) 
  
  seurat_merged <- IntegrateLayers(object = seurat_merged,
                                   method = method,
                                   orig.reduction = reduction,
                                   new.reduction = integration_tag,
                                   verbose = FALSE)
  
  seurat_merged <- seurat_clusters(seurat_merged,
                                  dim_n = dim_n,
                                  reduction = integration_tag,
                                  resolution = resolution,
                                  cluster_name = cluster_name)

  combinded_type <- seurat_merged@meta.data$condition
  combinded_type[grep("WT",combinded_type)] <- "WT"
  combinded_type[grep("KO",combinded_type)] <- "KO"
  combinded_type[grep("Het",combinded_type)] <- "Het"
  seurat_merged <- AddMetaData(seurat_merged,
                     metadata = combinded_type,
                     col.name = "condition_aggregated")
  return(seurat_merged)
}

#' Transfer Labels from Reference to Query
#' 
#' Uses Seurat's anchor-based transfer method to predict cell types in 
#' a query object using a reference object.
#' 
#' @param scrna The query Seurat object.
#' @param ref The reference Seurat object (must have 'cell_type' metadata).
#' @param dims Dimensions to use for anchor finding and data transfer.
#' @param reduction Reference reduction to use (default: "pca").
#'
#' @return The query Seurat object with predicted metadata.
#' @export
transfer_labels <- function(scrna,
                            ref,
                            dims = 1:30,
                            reduction = "pca") {
  anchors <- FindTransferAnchors(reference = ref,
                                 query = scrna,
                                 dims = dims,
                                 reference.reduction = reduction)
  
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
                      min_feature = min_features,
                      max_feature = max_features,
                      percent_mt = percent_mt,
                      nfeatures = n_var_features,
                      npcs = number_pcs)

seurat_list <- lapply(seurat_list,
                      FUN = seurat_clusters,
                      dim_n = number_pcs,
                      reduction = "pca", # PCA used as default
                      resolution = cluster_resolution)
saveRDS(seurat_list, file = "GAND_preprocessed.rds")

#-----------------------------------------------------------------------------#
# Integration
#-----------------------------------------------------------------------------#
seurat_integrated <- integrate_list(seurat_list,
                                    method = integration_method,
                                    reduction = "pca",
                                    integration_tag = integration_tag,
                                    dim_n = number_pcs,
                                    resolution = cluster_resolution,
                                    cluster_name = paste0(integration_tag,"_cluster_"))

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

