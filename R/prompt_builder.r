# prompt_builder.R
# Construct structured prompts for Ollama from DESeq2 and enrichment results

library(dplyr)
library(glue)

#' Format top DE genes into readable text block
#' @param sig Data frame. Significant DESeq2 results
#' @param top_n Integer. Number of genes to include
#' @param direction Character. "up", "down", or "both"
#' @return Character string
format_gene_block <- function(sig, top_n = 15, direction = "both") {
  
  df <- sig |>
    arrange(padj)
  
  if (direction == "up") {
    df <- df |> filter(log2FoldChange > 0)
  } else if (direction == "down") {
    df <- df |> filter(log2FoldChange < 0)
  }
  
  df <- df |> slice_head(n = top_n)
  
  if (nrow(df) == 0) return("None")
  
  df |>
    mutate(
      line = glue("{gene} (log2FC={round(log2FoldChange, 2)}, padj={signif(padj, 3)})")
    ) |>
    pull(line) |>
    paste(collapse = "\n")
}

#' Format pathway terms into readable text block
#' @param terms Character vector of pathway/GO terms
#' @param label Character. Section label
#' @return Character string
format_pathway_block <- function(terms, label = "Enriched pathways") {
  if (length(terms) == 0) return(glue("{label}: None detected"))
  glue("{label}:\n{paste(seq_along(terms), terms, sep='. ', collapse='\n')}")
}

#' Build summary statistics block
#' @param deseq_res List. Output from run_deseq2()
#' @return Character string
format_summary_block <- function(deseq_res) {
  sig   <- deseq_res$significant
  total <- nrow(deseq_res$results)
  n_sig <- nrow(sig)
  n_up  <- sum(sig$log2FoldChange > 0)
  n_dn  <- sum(sig$log2FoldChange < 0)
  
  glue(
    "Total genes tested: {total}
Significant DE genes: {n_sig} (padj < {deseq_res$padj_cutoff}, |log2FC| > {deseq_res$lfc_cutoff})
Upregulated: {n_up}
Downregulated: {n_dn}
Contrast: {deseq_res$contrast_level} vs {deseq_res$ref_level} (reference)"
  )
}

#' Build main interpretation prompt
#' @param deseq_res List. Output from run_deseq2()
#' @param pathway_terms List. Output from extract_top_terms()
#' @param top_n Integer. Number of top genes to include
#' @return Character string
build_interpretation_prompt <- function(deseq_res,
                                         pathway_terms,
                                         top_n = 15) {
  
  summary_block  <- format_summary_block(deseq_res)
  up_block       <- format_gene_block(deseq_res$significant, top_n, "up")
  dn_block       <- format_gene_block(deseq_res$significant, top_n, "down")
  go_block       <- format_pathway_block(pathway_terms$go_terms,   "GO Biological Process terms")
  kegg_block     <- format_pathway_block(pathway_terms$kegg_terms, "KEGG pathways")
  
  glue(
    "You are an expert bioinformatics analyst. Interpret the following differential 
expression results for a human RNA-seq experiment.

--- EXPERIMENT SUMMARY ---
{summary_block}

--- TOP UPREGULATED GENES ---
{up_block}

--- TOP DOWNREGULATED GENES ---
{dn_block}

--- PATHWAY ENRICHMENT ---
{go_block}

{kegg_block}

--- INSTRUCTIONS ---
1. Summarise the key biological processes altered between {deseq_res$contrast_level} 
   and {deseq_res$ref_level} in 3-5 sentences.
2. Highlight the most biologically significant gene-pathway associations.
3. Note any potentially interesting or unexpected findings.
4. Do not speculate beyond the data provided.
5. Do not fabricate gene functions or pathway associations.
6. Use precise scientific language appropriate for a research audience.
7. Structure your response with these sections:
   - Biological Summary
   - Key Gene-Pathway Associations
   - Notable Findings"
  )
}

#' Build upregulated genes specific prompt
#' @param deseq_res List. Output from run_deseq2()
#' @param top_n Integer
#' @return Character string
build_upregulated_prompt <- function(deseq_res, top_n = 15) {
  up_block <- format_gene_block(deseq_res$significant, top_n, "up")
  
  glue(
    "The following genes are significantly upregulated in {deseq_res$contrast_level} 
vs {deseq_res$ref_level} in a human RNA-seq experiment:

{up_block}

Provide a concise biological interpretation of these upregulated genes in 2-3 sentences.
Focus on shared biological functions or pathways.
Do not speculate beyond established gene functions.
Do not fabricate associations."
  )
}

#' Build downregulated genes specific prompt
#' @param deseq_res List. Output from run_deseq2()
#' @param top_n Integer
#' @return Character string
build_downregulated_prompt <- function(deseq_res, top_n = 15) {
  dn_block <- format_gene_block(deseq_res$significant, top_n, "down")
  
  glue(
    "The following genes are significantly downregulated in {deseq_res$contrast_level} 
vs {deseq_res$ref_level} in a human RNA-seq experiment:

{dn_block}

Provide a concise biological interpretation of these downregulated genes in 2-3 sentences.
Focus on shared biological functions or pathways.
Do not speculate beyond established gene functions.
Do not fabricate associations."
  )
}
