#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(liana)
  library(readr)
  library(dplyr)
  library(tibble)
})

counts_file      <- "counts_OldvsYoung.csv"
meta_file        <- "meta_OldvsYoung.csv"
cluster_map_file <- "Old_cluster_map.csv"
output_file      <- "liana_results_old.csv"
normalize_input  <- TRUE

read_counts_as_matrix <- function(path) {
  counts_df <- read_csv(path, show_col_types = FALSE)
  if (!"Barcodes" %in% colnames(counts_df)) {
    stop("The counts file must contain a 'Barcodes' column.")
  }
  barcodes <- counts_df$Barcodes
  counts_df <- counts_df %>% select(-Barcodes)
  counts_mat <- as.matrix(counts_df)
  mode(counts_mat) <- "numeric"
  rownames(counts_mat) <- barcodes
  counts_mat <- t(counts_mat)
  return(counts_mat)
}

prepare_seurat_object <- function(counts_file, meta_file, cluster_map_file) {
  counts_mat <- read_counts_as_matrix(counts_file)
  meta_df <- read_csv(meta_file, show_col_types = FALSE)
  cluster_map <- read_csv(cluster_map_file, show_col_types = FALSE)

  required_meta_cols <- c("Barcodes", "Cluster", "source")
  if (!all(required_meta_cols %in% colnames(meta_df))) {
    stop("Metadata file must contain columns: ", paste(required_meta_cols, collapse = ", "))
  }
  required_map_cols <- c("cluster_id", "cell_type")
  if (!all(required_map_cols %in% colnames(cluster_map))) {
    stop("Cluster map file must contain columns: ", paste(required_map_cols, collapse = ", "))
  }

  meta_df <- meta_df %>% left_join(cluster_map, by = c("Cluster" = "cluster_id"))

  if (any(is.na(meta_df$cell_type))) {
    unmapped <- unique(meta_df$Cluster[is.na(meta_df$cell_type)])
    stop("Some clusters were not mapped to cell types: ", paste(unmapped, collapse = ", "))
  }

  meta_df <- as.data.frame(meta_df)
  rownames(meta_df) <- meta_df$Barcodes

  shared_cells <- intersect(colnames(counts_mat), rownames(meta_df))
  if (length(shared_cells) == 0) {
    stop("No overlapping barcodes were found between counts and metadata.")
  }

  counts_mat <- counts_mat[, shared_cells, drop = FALSE]
  meta_df <- meta_df[shared_cells, , drop = FALSE]

  seu <- CreateSeuratObject(counts = counts_mat, meta.data = meta_df)

  if (normalize_input) {
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  }

  return(seu)
}

run_liana_and_save <- function(seu_obj, output_file) {
  Idents(seu_obj) <- seu_obj@meta.data$cell_type

  res_list <- liana_wrap(
    seu_obj,
    method = c("cellphonedb", "connectome", "natmi", "sca", "logfc"),
    resource = "Consensus",
    expr_prop = 0.1,
    min_cells = 5,
    permutation.params = list(nperms = 1000),
    assay = DefaultAssay(seu_obj)
  )

  res_agg <- liana_aggregate(res_list)

  res_tbl <- as_tibble(res_agg) %>%
    select(
      source, target, ligand_complex, receptor_complex,
      lr_means, cellphone_pvals, expr_prod, scaled_weight,
      lr_logfc, spec_weight, lrscore, specificity_rank, magnitude_rank
    )

  write_csv(res_tbl, output_file)
  message("Saved LIANA result table to: ", output_file)
}

seu <- prepare_seurat_object(counts_file, meta_file, cluster_map_file)
old_obj <- subset(seu, subset = source == "Old")
if (ncol(old_obj) == 0) {
  stop("No cells found for source == 'Old'.")
}
run_liana_and_save(old_obj, output_file)
