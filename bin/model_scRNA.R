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
# Parse arguments
argv <- parse_args(p)

print(argv)

# Assign to variables
annotated <- argv$annotated
gene_sets <- fromJSON(argv$gene_sets)
mut_genes <- fromJSON(argv$mut_genes)
score_type <- fromJSON(argv$score_type)


max_size <- 100000 * 1024^2
options(future.globals.maxSize = max_size)
#-----------------------------------------------------------------------------#
# UTILS
#-----------------------------------------------------------------------------#
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
  #table(meta_data$sample, meta_data$condition_aggregated)
  cell_counts <- table(meta_data$cell_type, meta_data$sample)
  dropped_cell <- rownames(cell_counts)[apply(cell_counts,1,function(x,mc){return(any(x < mc))},mc = min_cells)]
  # remove - not enough cells to get any statistical insights
  meta_data <- meta_data %>% filter(!cell_type %in% dropped_cell)


  model_main <- lmer(
    gene_set1 ~ condition_aggregated * cell_type + (1 | sample),
    data = meta_data,
    REML = TRUE
  )

  # Get p-values and confidence intervals using emmeans (estimated marginal means)
  # Test 1: Overall condition effect (averaged across cell types)
  condition_effect <- emmeans(model_main, ~ condition_aggregated)
  

  # Test 2: Condition effect within each cell type (most relevant for your question)
  condition_by_celltype <- emmeans(model_main, ~ condition_aggregated | cell_type)


  # Pairwise comparisons: WT vs Condition1 vs Condition2 within each cell type
  pairwise_results <- pairs(condition_by_celltype, adjust = "BH")

  # Extract all p-values from pairwise comparisons
  pairwise_df <- as.data.frame(pairwise_results)
  pairwise_df$p_value_adjusted <- p.adjust(pairwise_df$p.value, method = "BH")
  pairwise_df$significant <- pairwise_df$p_value_adjusted < 0.05

  # Preparing plotting table.
  emmeans_df <- as.data.frame(condition_by_celltype)

  # Prepare pairwise significance annotations (WT vs others)
  pairwise_df <- as.data.frame(pairwise_results) %>%
    filter(grepl("WT -", contrast)) %>%
    mutate(
      group1 = gsub(" - .*", "", contrast),
      group2 = gsub(".* - ", "", contrast),
      # Assign stars based on p-value
      sig_label = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    ) %>%
    select(cell_type, condition_aggregated = group2, sig_label)

  # Merge with emmeans
  plot_df <- emmeans_df %>%
    left_join(pairwise_df, by = c("cell_type", "condition_aggregated"))
  # Merge with gene set names
  gene_set <- paste0(gene_set,collapse = "|")
  plot_df$gene_set <- gene_set
  return(plot_df)
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
    write.csv(tmp, file = paste0(score_type[s],"_",g,"_geneset_list.csv"))
  }
  
}

#-----------------------------------------------------------------------------#
# DONE
#-----------------------------------------------------------------------------#

