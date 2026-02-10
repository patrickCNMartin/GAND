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
p <- arg_parser("Model GAND Cells")

# Add arguments
p <- add_argument(p, "--annotated", short = "-a", help = "Annotated seurat objecy", type = "character")
p <- add_argument(p, "--gene_sets", short = "-g",help = "Gene sets to check for enrichement", type = "character")
p <- add_argument(p, "--mut_genes", short = "-m",help = "Mutally Excuslive Genes", type = "character")
p <- add_argument(p, "--score_type", short = "-s", help = "How to score gene sets for LLM modelling",type = "character")
p <- add_argument(p, "--min_cells", short = "-mc", help = "If cell type has less then min number of cells - removed from modelling",type = "numeric")
# Parse arguments
argv <- parse_args(p)

print(argv)

# Assign to variables
annotated <- argv$annotated
gene_sets <- fromJSON(argv$gene_sets)
mut_genes <- fromJSON(argv$mut_genes)
score_type <- fromJSON(argv$score_type)
min_cells <- argv$min_cells


max_size <- 100000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# UTILS
#-----------------------------------------------------------------------------#
#' Mutual Expression of Two Genes
#'
#' Extracts expression data for two specific genes across all cells and merges 
#' it with key metadata (cell type, sample, and condition) for downstream 
#' mutual exclusivity or co-expression analysis.
#'
#' @param seurat_annotated A Seurat object containing the gene expression data and metadata.
#' @param genes A character vector of exactly two gene names to extract.
#'
#' @return A data frame where each row is a cell, containing gene expression values 
#'         and columns for cell_type, sample, and condition_aggregated.
#' @export
#'
#' @import Seurat
mut_ex_genes <- function(seurat_annotated,
                         genes){
  if (length(genes) > 2) {
    stop("Can only check 2 genes at a time.")
  }
  seurat_list <- SplitObject(seurat_annotated, split.by = "sample")
  counts <- lapply(seurat_list, GetAssayData, layer = "data")
  counts <- do.call("cbind", counts)
  ex_genes <- t(counts[genes,])
  ex_genes <- data.frame("cell_id" = rownames(ex_genes),as.matrix(ex_genes))

  cell_type <- seurat_annotated@meta.data$predicted.id
  names(cell_type) <- rownames(seurat_annotated@meta.data)
  cell_type <- cell_type[match(ex_genes$cell_id,names(cell_type))]

  sample <- seurat_annotated@meta.data$sample
  names(sample) <- rownames(seurat_annotated@meta.data)
  sample <- sample[match(ex_genes$cell_id,names(sample))]

  condition_aggregated <- seurat_annotated@meta.data$condition_aggregated
  names(condition_aggregated) <- rownames(seurat_annotated@meta.data)
  condition_aggregated <- condition_aggregated[match(ex_genes$cell_id,names(condition_aggregated))]

  ex_genes$cell_type <- cell_type
  ex_genes$sample <- sample
  ex_genes$condition_aggregated <- condition_aggregated
  return(ex_genes)
}

#' Add Aggregate Count Score to Metadata
#'
#' Calculates the sum of expression for a set of features per cell across 
#' samples and adds the result to the Seurat object's metadata.
#'
#' @param seurat_annotated A Seurat object.
#' @param features A character vector of gene names to sum.
#' @param name Character string used as the prefix for the new metadata column.
#'
#' @return The Seurat object with an additional metadata column named `paste0(name, "1")`.
#' @export
#'
#' @import Seurat
AddCountScore <- function(seurat_annotated,
                          features,
                          name) {
  seurat_list <- SplitObject(seurat_annotated, split.by = "sample")
  counts <- lapply(seurat_list, GetAssayData, layer = "data")
  counts <- lapply(counts, function(s, gs,name){
    sub_counts <- s[gs,]
    if (!is.null(dim(sub_counts))){
      sum_counts <- as.vector(apply(sub_counts,2,sum))
      df <- data.frame(colnames(sub_counts),sum_counts)
    } else {
      sum_counts <- as.vector(sub_counts)
      df <- data.frame(names(sub_counts),sum_counts)
    }
    colnames(df) <- c("cell_id", paste0(name,"1"))
    return(df)
  }, gs = features, name = name)
  counts <- do.call("rbind",counts)
  row_order <- rownames(seurat_annotated@meta.data)
  counts <- as.vector(counts[match(row_order,counts$cell_id),paste0(name,"1")])
  seurat_annotated <- AddMetaData(seurat_annotated, counts, col.name = paste0(name,"1"))
  return(seurat_annotated)
}
#' Gene Set Enrichment Analysis using Mixed Effects Models
#'
#' This function calculates gene set scores for single-cell data and performs 
#' statistical testing using a linear mixed-effects model to account for 
#' sample-level variation.
#'
#' @param seurat_annotated A Seurat object with metadata for 'sample', 'condition_aggregated', and 'predicted.id'.
#' @param gene_set A character vector or list of genes to score.
#' @param score_type Character, either "module" (Seurat::AddModuleScore) or "counts" (AddCountScore).
#' @param seed Numeric, seed for reproducibility in module scoring.
#' @param min_cells Integer, minimum cells per sample for a cell type to be included in the model.
#'
#' @return A data frame containing estimated marginal means (emmeans), standard errors, 
#'         p-values for comparisons against WT, and significance labels.
#' @export
#'
#' @import Seurat
#' @import lme4
#' @import emmeans
#' @import dplyr
gene_set_enrichement <- function(seurat_annotated,
                                 gene_set,
                                 score_type = "module",
                                 seed = 42,
                                 min_cells = 1){

  # Add module score using Seurat
  if (score_type == "module") {
    seurat_obj <- AddModuleScore(seurat_annotated,
                               features = gene_set,
                               name = "gene_set",
                               seed = seed)
  } else if (score_type == "counts") {
    seurat_obj <- AddCountScore(seurat_annotated,
                               features = gene_set,
                               name = "gene_set")
  } else {
    stop("Unknown score type request")
  }
  
  # Extract metadata with scores
  meta_data <- seurat_obj@meta.data[,c("sample","condition","condition_aggregated","predicted.id", "gene_set1")]
  meta_data$cell_barcode <- rownames(meta_data)
  colnames(meta_data) <- gsub("predicted.id","cell_type",colnames(meta_data))

  meta_data$sample <- factor(meta_data$sample)
  meta_data$condition_aggregated <- factor(meta_data$condition_aggregated, 
                               levels = c("WT", "Het", "KO"))
  meta_data$cell_type <- factor(meta_data$cell_type)

  # Filter out cell types with insufficient cell counts per sample
  cell_counts <- table(meta_data$cell_type, meta_data$sample)
  dropped_cell <- rownames(cell_counts)[apply(cell_counts, 1, function(x, mc){return(any(x < mc))}, mc = min_cells)]
  meta_data <- meta_data %>% filter(!cell_type %in% dropped_cell)

  # Mixed Effects Model
  model_main <- lmer(
    gene_set1 ~ condition_aggregated * cell_type + (1 | sample),
    data = meta_data,
    REML = TRUE
  )

  # Condition effect within each cell type
  condition_by_celltype <- emmeans(model_main, ~ condition_aggregated | cell_type)

  # Pairwise comparisons: WT vs Others within each cell type
  pairwise_results <- pairs(condition_by_celltype, adjust = "BH")
  
  # Process results into a clean dataframe
  pairwise_df <- as.data.frame(pairwise_results) %>%
    filter(grepl("WT -", contrast)) %>%
    mutate(
      group1 = gsub(" - .*", "", contrast),
      group2 = gsub(".* - ", "", contrast),
      sig_label = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    ) %>%
    select(cell_type, condition_aggregated = group2, p.value, sig_label)

  # Merge means and significance
  emmeans_df <- as.data.frame(condition_by_celltype)
  results_df <- emmeans_df %>%
    left_join(pairwise_df, by = c("cell_type", "condition_aggregated"))

  # Attach gene set metadata
  results_df$gene_set <- paste0(unlist(gene_set), collapse = "|")
  
  return(results_df)
}

#-----------------------------------------------------------------------------#
# DATA LOADING
#-----------------------------------------------------------------------------#
annotated <- readRDS(annotated)
#-----------------------------------------------------------------------------#
# Fetch Annotations
#-----------------------------------------------------------------------------#
df <- FetchData(annotated,
                c("umap_1",
                  "umap_2",
                  "predicted.id",
                  "sample",
                  "condition",
                  "condition_aggregated"))
#save plotting df as source data
write.csv(df, file = "GAND_seurat_annotated.csv")
#-----------------------------------------------------------------------------#
# Mutal exclusive genes
#-----------------------------------------------------------------------------#
ex_genes <- mut_ex_genes(annotated,
                         mut_genes)

write.csv(ex_genes, file = "mutually_exclusive_genes.csv")
#-----------------------------------------------------------------------------#
# Gene set enrichement by cell and condition
#-----------------------------------------------------------------------------#
for (g in seq_along(gene_sets)) {
  for (s in seq_along(score_type)) {
    tmp <- gene_set_enrichement(annotated,
                                gene_sets[[g]],
                                score_type = score_type[s],
                                seed = 42,
                                min_cells = 1)
    write.csv(tmp, file = paste0(score_type[s], "_", g, "_geneset_list.csv"))
  }
  
}

#-----------------------------------------------------------------------------#
# DONE
#-----------------------------------------------------------------------------#

