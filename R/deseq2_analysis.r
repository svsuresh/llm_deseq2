# deseq2_analysis.R
# DESeq2 wrapper functions for Shiny app

library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)


#' Validate counts matrix
#' @param counts Data frame. Raw counts matrix
#' @param metadata Data frame. Sample metadata
#' @param condition Character. Condition column name
#' @return List with valid flag and message
validate_inputs <- function(counts, metadata, condition) {
  
  # Check condition column exists
  if (!condition %in% colnames(metadata)) {
    return(list(valid = FALSE, 
                msg = paste("Condition column", condition, "not found in metadata")))
  }
  
  # Check sample names match
  meta_samples <- metadata[[1]]
  count_samples <- colnames(counts)[-1]
  
  if (!all(meta_samples %in% count_samples)) {
    missing <- meta_samples[!meta_samples %in% count_samples]
    return(list(valid = FALSE,
                msg = paste("Samples in metadata not found in counts:", 
                            paste(missing, collapse = ", "))))
  }
  
  # Check exactly two levels
  levels <- unique(metadata[[condition]])
  if (length(levels) != 2) {
    return(list(valid = FALSE,
                msg = paste("Condition must have exactly 2 levels. Found:", 
                            paste(levels, collapse = ", "))))
  }
  
  # Check minimum samples per group
  counts_per_group <- table(metadata[[condition]])
  if (any(counts_per_group < 2)) {
    return(list(valid = FALSE,
                msg = "Each condition group must have at least 2 samples"))
  }
  
  list(valid = TRUE, msg = "OK")
}

#' Prepare DESeq2 dataset
#' @param counts Data frame. Raw counts with gene column first
#' @param metadata Data frame. Sample metadata
#' @param condition Character. Condition column name
#' @param ref_level Character. Reference level for contrast
#' @return DESeqDataSet object
prepare_dds <- function(counts, metadata, condition, ref_level) {
  
  # Format counts matrix
  count_mat <- counts |>
    column_to_rownames(colnames(counts)[1]) |>
    as.matrix()
  
  # Align samples
  meta_df <- metadata |>
    column_to_rownames(colnames(metadata)[1]) |>
    as.data.frame()
  
  count_mat <- count_mat[, rownames(meta_df)]
  
  # Set reference level
  meta_df[[condition]] <- relevel(
    factor(meta_df[[condition]]), 
    ref = ref_level
  )
  
  # Build DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = meta_df,
    design    = as.formula(paste("~", condition))
  )
  
  # Pre-filter low count genes
  dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
  
  dds
}

#' Run DESeq2 analysis
#' @param dds DESeqDataSet object
#' @param condition Character. Condition column name
#' @param ref_level Character. Reference level
#' @param padj_cutoff Numeric. Adjusted p-value cutoff
#' @param lfc_cutoff Numeric. Log2 fold change cutoff
#' @return List with dds, results, and significant genes
run_deseq2 <- function(dds, 
                        condition, 
                        ref_level,
                        padj_cutoff = 0.05, 
                        lfc_cutoff  = 1) {
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get contrast levels
  lvls <- levels(dds[[condition]])
  contrast_level <- lvls[lvls != ref_level]
  
  # Extract results
  res <- results(
    dds,
    contrast    = c(condition, contrast_level, ref_level),
    alpha       = padj_cutoff
  ) |>
    as.data.frame() |>
    rownames_to_column("gene") |>
    arrange(padj) |>
    filter(!is.na(padj), !is.na(log2FoldChange))
  
  # Significant genes
  sig <- res |>
    filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
  
  list(
    dds            = dds,
    results        = res,
    significant    = sig,
    contrast_level = contrast_level,
    ref_level      = ref_level,
    padj_cutoff    = padj_cutoff,
    lfc_cutoff     = lfc_cutoff
  )
}

#' Generate volcano plot
#' @param res Data frame. Full DESeq2 results
#' @param padj_cutoff Numeric. Adjusted p-value cutoff for the volcano plot
#' @param lfc_up Numeric. Positive log2FC cutoff (must be > 0)
#' @param lfc_dn Numeric. Negative log2FC cutoff (must be < 0)
#' @param top_n Integer. Number of top genes to label
#' @return ggplot object
plot_volcano <- function(res,
                         padj_cutoff = 0.05,
                         lfc_up      = 1,
                         lfc_dn      = -1,
                         top_n       = 15) {

  res <- res |>
    mutate(
      neg_log10_padj = -log10(padj),
      status = case_when(
        padj < padj_cutoff & log2FoldChange >  lfc_up ~ "Up",
        padj < padj_cutoff & log2FoldChange <  lfc_dn ~ "Down",
        TRUE ~ "NS"
      )
    )

  # Top genes to label
  top_genes <- res |>
    filter(status != "NS") |>
    slice_head(n = top_n)

  ggplot(res, aes(log2FoldChange, neg_log10_padj, color = status)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_text_repel(
      data    = top_genes,
      aes(label = gene),
      size    = 3,
      max.overlaps = 20
    ) +
    geom_vline(xintercept = c(lfc_dn, lfc_up),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(Up = "#E64B35", Down = "#4DBBD5", NS = "grey70")) +
    labs(
      x     = "Log2 Fold Change",
      y     = "-Log10 Adjusted P-value",
      color = "Status",
      title = "Volcano Plot"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "top",
      text            = element_text(color = "black"),
      axis.text       = element_text(color = "black")
    )
}

#' Generate PCA plot
#' @param dds DESeqDataSet object
#' @param condition Character. Condition column name
#' @return ggplot object
plot_pca <- function(dds, condition) {
  vsd <- vst(dds, blind = TRUE)
  plotPCA(vsd, intgroup = condition) +
    theme_bw(base_size = 13) +
    ggtitle("PCA Plot") +
    theme(
      legend.position = "top",
      text            = element_text(color = "black"),
      axis.text       = element_text(color = "black")
    )
}
