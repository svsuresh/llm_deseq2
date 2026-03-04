# enrichment.R
# Pathway enrichment analysis using clusterProfiler (Human only)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(ggplot2)


#' Convert gene symbols to Entrez IDs
#' @param genes Character vector of gene symbols
#' @return Named character vector of Entrez IDs
symbols_to_entrez <- function(genes) {
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys    = genes,
    column  = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  # Remove NAs
  mapped <- mapped[!is.na(mapped)]
  as.character(mapped)
}

#' Run GO biological process enrichment
#' @param sig_genes Character vector. Significant gene symbols
#' @param all_genes Character vector. Background gene symbols
#' @param padj_cutoff Numeric. Adjusted p-value cutoff
#' @param min_gs Integer. Minimum gene set size
#' @param max_gs Integer. Maximum gene set size
#' @return enrichResult object or NULL
run_go_enrichment <- function(sig_genes,
                               all_genes,
                               padj_cutoff = 0.05,
                               min_gs      = 10,
                               max_gs      = 500) {
  
  sig_entrez <- symbols_to_entrez(sig_genes)
  bg_entrez  <- symbols_to_entrez(all_genes)
  
  if (length(sig_entrez) < 5) {
    return(NULL)
  }
  
  ego <- enrichGO(
    gene          = sig_entrez,
    universe      = bg_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = padj_cutoff,
    qvalueCutoff  = 0.2,
    minGSSize     = min_gs,
    maxGSSize     = max_gs,
    readable      = TRUE
  )
  
  if (is.null(ego) || nrow(ego@result) == 0) {
    return(NULL)
  }
  
  ego
}

#' Run KEGG pathway enrichment
#' @param sig_genes Character vector. Significant gene symbols
#' @param all_genes Character vector. Background gene symbols
#' @param padj_cutoff Numeric
#' @return enrichResult object or NULL
run_kegg_enrichment <- function(sig_genes,
                                 all_genes,
                                 padj_cutoff = 0.05) {
  
  sig_entrez <- symbols_to_entrez(sig_genes)
  bg_entrez  <- symbols_to_entrez(all_genes)
  
  if (length(sig_entrez) < 5) {
    return(NULL)
  }
  
  ekegg <- enrichKEGG(
    gene          = sig_entrez,
    universe      = bg_entrez,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = padj_cutoff
  )
  
  if (is.null(ekegg) || nrow(ekegg@result) == 0) {
    return(NULL)
  }
  
  ekegg
}

#' Extract top pathway terms for LLM prompt
#' @param ego enrichResult object. GO results
#' @param ekegg enrichResult object. KEGG results
#' @param top_n Integer. Number of top terms to extract
#' @return List with go_terms and kegg_terms character vectors
extract_top_terms <- function(ego   = NULL,
                               ekegg = NULL,
                               top_n = 10) {
  
  go_terms <- character(0)
  if (!is.null(ego)) {
    go_terms <- ego@result |>
      filter(p.adjust < 0.05) |>
      arrange(p.adjust) |>
      slice_head(n = top_n) |>
      pull(Description)
  }
  
  kegg_terms <- character(0)
  if (!is.null(ekegg)) {
    kegg_terms <- ekegg@result |>
      filter(p.adjust < 0.05) |>
      arrange(p.adjust) |>
      slice_head(n = top_n) |>
      pull(Description)
  }
  
  list(go_terms = go_terms, kegg_terms = kegg_terms)
}

#' Plot GO dotplot
#' @param ego enrichResult object
#' @param top_n Integer. Number of terms to show
#' @return ggplot object or NULL
plot_go_dotplot <- function(ego, top_n = 20) {
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  
  dotplot(ego, showCategory = top_n) +
    ggtitle("GO Biological Process Enrichment") +
    theme_bw(base_size = 12) +
    theme(
      text        = element_text(color = "black"),
      axis.text   = element_text(color = "black"),
      axis.text.y = element_text(size = 9, color = "black")
    )
}

#' Plot KEGG dotplot
#' @param ekegg enrichResult object
#' @param top_n Integer. Number of terms to show
#' @return ggplot object or NULL
plot_kegg_dotplot <- function(ekegg, top_n = 20) {
  if (is.null(ekegg) || nrow(ekegg@result) == 0) return(NULL)

  dotplot(ekegg, showCategory = top_n) +
    ggtitle("KEGG Pathway Enrichment") +
    theme_bw(base_size = 12) +
    theme(
      text        = element_text(color = "black"),
      axis.text   = element_text(color = "black"),
      axis.text.y = element_text(size = 9, color = "black")
    )
}

#' Plot GO emap (enrichment map)
#' @param ego enrichResult object
#' @param top_n Integer
#' @return ggplot object or NULL
plot_go_emap <- function(ego, top_n = 30) {
  if (is.null(ego) || nrow(ego@result) < 2) return(NULL)

  tryCatch({
    ego2 <- pairwise_termsim(ego)
    emapplot(ego2, showCategory = top_n) +
      ggtitle("GO Term Enrichment Map") +
      theme(
        text      = element_text(color = "black"),
        axis.text = element_text(color = "black")
      )
  }, error = function(e) NULL)
}
