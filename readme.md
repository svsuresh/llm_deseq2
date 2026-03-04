# RNA-seq LLM Interpreter

A Shiny application for differential expression analysis of human RNA-seq data with locally deployed LLM-assisted biological interpretation.

Combines DESeq2, clusterProfiler, and a locally running Ollama model to produce an interpretable HTML or PDF report — without sending any data to external APIs.

---

## What It Does

1. Accepts user-supplied raw counts and metadata files
2. Runs DESeq2 differential expression analysis
3. Runs GO and KEGG pathway enrichment
4. Queries a locally deployed LLM (via Ollama) for biological interpretation
5. Renders a self-contained HTML or PDF report

All analysis runs locally. No data leaves your machine.

---

## Requirements

### System
- R >= 4.3
- [Ollama](https://ollama.com) installed and running
- [Quarto CLI](https://quarto.org/docs/get-started/) installed

### R Packages

```r
# CRAN
install.packages(c(
  "shiny", "bslib", "httr2", "dplyr", "tibble",
  "ggplot2", "ggrepel", "glue", "here", "quarto"
))

# Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
))
```

### Ollama Setup

```bash
# Start Ollama
ollama serve

# Pull a model (choose one)
ollama pull mistral        # recommended, good balance of speed and quality
ollama pull llama3.2       # lighter, faster
ollama pull llama3.1:8b    # more capable, slower
```

---

## Installation

```bash
git clone https://github.com/svsuresh/rnaseq-llm-interpreter.git
cd rnaseq-llm-interpreter
```

---

## Running the App

```r
shiny::runApp(".")
```

Or from terminal:

```bash
Rscript -e "shiny::runApp('.')"
```

---

## Input File Format

Both files must be CSV or TSV with a header row.

### Counts Matrix
- First column: gene symbols (HGNC)
- Remaining columns: one per sample
- Values: raw integer counts (not normalised, not log-transformed)

```
gene,sample1,sample2,sample3,sample4
BRCA1,100,200,150,300
TP53,50,80,60,90
```

### Metadata
- First column: sample names matching counts column headers exactly
- Remaining columns: sample attributes
- Condition column must have exactly 2 levels

```
sample,condition
sample1,control
sample2,control
sample3,treated
sample4,treated
```

Example files are provided in `example_data/`.

---

## Usage

1. Open the app — check the **Setup** tab to confirm Ollama is detected
2. Upload counts matrix and metadata file
3. Select the condition column and reference level
4. Adjust cutoffs if needed (default: padj < 0.05, |log2FC| > 1)
5. Select your Ollama model
6. Click **Run Analysis**
7. Monitor progress in the **Log** tab
8. Download the report (HTML or PDF) when analysis completes

---

## Output Report Sections

- Experiment summary
- PCA plot (quality control)
- Volcano plot
- Top upregulated and downregulated gene tables
- Full results table
- GO Biological Process enrichment (dotplot + enrichment map)
- KEGG pathway enrichment
- LLM-generated biological interpretation (overall, upregulated, downregulated)
- Session information

---

## A Note on LLM Interpretation

LLM-generated sections are clearly labelled in the report. Outputs are based strictly on the differential expression and pathway enrichment results — the model is explicitly instructed not to speculate beyond the data provided.

All LLM outputs in this application were reviewed for biological accuracy before inclusion in the report. Users are advised to critically evaluate interpretations in the context of their experimental system. LLM outputs are a starting point for interpretation, not a substitute for domain expertise.

---

## Limitations

- Human samples only (hg38, HGNC gene symbols)
- Single condition with two levels
- Requires Ollama running locally on port 11434
- Report rendering requires Quarto CLI

---

## Troubleshooting

| Problem | Fix |
|--------|-----|
| Ollama not detected | Run `ollama serve` in terminal before launching app |
| No models in dropdown | Run `ollama pull mistral` |
| Gene ID conversion fails | Ensure gene column contains HGNC symbols, not Ensembl IDs |
| Report render fails | Verify Quarto CLI is installed and on PATH |
| Ollama timeout | Switch to a smaller/faster model like `llama3.2` |
| DESeq2 convergence warning | Normal for some datasets, results are still valid |

---

## Repository Structure

```
rnaseq-llm-interpreter/
├── app.R                  # Shiny application
├── R/
│   ├── ollama.R           # Ollama API bridge
│   ├── deseq2_analysis.R  # DESeq2 wrapper functions
│   ├── enrichment.R       # GO and KEGG enrichment
│   └── prompt_builder.R   # LLM prompt construction
├── report/
│   └── report.qmd         # Quarto report template
├── example_data/
│   ├── counts.csv         # Example counts matrix
│   └── metadata.csv       # Example metadata
└── README.md
```

---

## License

MIT

---

## Author

Venkata Suresh Kumar, PhD  
[github.com/svsuresh](https://github.com/svsuresh) · [linkedin.com/in/svsk](https://linkedin.com/in/svsk)
