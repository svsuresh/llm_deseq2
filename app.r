# app.R
# RNA-seq LLM Interpreter — Shiny Application

library(shiny)
library(bslib)
library(dplyr)
library(here)

source(here("R", "ollama.r"))
source(here("R", "deseq2_analysis.r"))
source(here("R", "enrichment.r"))
source(here("R", "prompt_builder.r"))

# -------------------------------------------------------------------------
# UI
# -------------------------------------------------------------------------

ui <- page_sidebar(
  title = "Transcriptomics LLM Interpreter",
  theme = bs_theme(bootswatch = "flatly", base_font = font_google("Inter")),

  sidebar = sidebar(
    width = 320,

    # --- Inputs ---
    h5("1. Upload Data", class = "mt-2 fw-bold"),

    fileInput(
      "counts_file", "Raw Counts Matrix",
      accept  = c(".csv", ".tsv", ".txt", ".gz"),
      placeholder = "CSV or TSV (plain or gzipped)"
    ),

    fileInput(
      "meta_file", "Sample Metadata",
      accept  = c(".csv", ".tsv", ".txt", ".gz"),
      placeholder = "CSV or TSV (plain or gzipped)"
    ),

    hr(),
    h5("2. Analysis Parameters", class = "fw-bold"),

    uiOutput("condition_select"),
    uiOutput("ref_level_select"),

    numericInput(
      "padj_cutoff", "Adj. P-value Cutoff",
      value = 0.05, min = 0.001, max = 0.1, step = 0.01
    ),

    numericInput(
      "lfc_cutoff", "Log2FC Cutoff",
      value = 1, min = 0.5, max = 3, step = 0.5
    ),

    hr(),
    h5("3. LLM Settings", class = "fw-bold"),

    uiOutput("model_select"),

    numericInput(
      "temperature", "Temperature",
      value = 0.3, min = 0.0, max = 1.0, step = 0.1
    ),

    hr(),

    actionButton(
      "run_btn", "Run Analysis",
      class = "btn-primary w-100",
      icon  = icon("play")
    ),

    br(), br(),

    uiOutput("download_ui")
  ),

  # --- Main panel ---
  navset_tab(

    nav_panel(
      "Setup",
      icon = icon("circle-info"),
      card(
        card_header("Getting Started"),
        card_body(
          h5("Requirements"),
          tags$ul(
            tags$li("Ollama must be running locally on port 11434"),
            tags$li("At least one model must be pulled (e.g. ", tags$code("ollama pull mistral"), ")"),
            tags$li("R packages: DESeq2, clusterProfiler, org.Hs.eg.db, enrichplot, httr2, quarto")
          ),
          h5("Input File Format"),
          h6("Counts Matrix (CSV/TSV)"),
          tags$pre(
"gene,sample1,sample2,sample3,sample4
BRCA1,100,200,150,300
TP53,50,80,60,90"
          ),
          h6("Metadata (CSV/TSV)"),
          tags$pre(
"sample,condition
sample1,control
sample2,control
sample3,treated
sample4,treated"
          ),
          h5("Notes"),
          tags$ul(
            tags$li("First column of counts must be gene symbols"),
            tags$li("First column of metadata must be sample names matching counts columns"),
            tags$li("Condition must have exactly 2 levels"),
            tags$li("Minimum 2 samples per group required"),
            tags$li("Analysis is for human samples only (hg38)"),
            tags$li("Gzipped CSV files (.csv.gz) are supported")
          )
        )
      )
    ),

    nav_panel(
      "QC",
      icon = icon("magnifying-glass-chart"),
      card(
        card_header("Principal Component Analysis"),
        card_body(plotOutput("pca_plot", height = "500px"))
      )
    ),

    nav_panel(
      "Differential Expression",
      icon = icon("chart-line"),
      card(
        card_header("Volcano Plot"),
        card_body(
          layout_columns(
            col_widths = c(4, 4, 4),
            numericInput(
              "volcano_lfc_up", "Upregulated LFC Cutoff",
              value = 1, min = 0, max = 10, step = 0.5
            ),
            numericInput(
              "volcano_lfc_dn", "Downregulated LFC Cutoff",
              value = -1, min = -10, max = 0, step = 0.5
            ),
            numericInput(
              "volcano_padj", "Volcano P-value Cutoff",
              value = 0.05, min = 0.0001, max = 0.1, step = 0.01
            )
          ),
          plotOutput("volcano_plot", height = "500px")
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Top Upregulated Genes"),
          card_body(tableOutput("up_table"))
        ),
        card(
          card_header("Top Downregulated Genes"),
          card_body(tableOutput("dn_table"))
        )
      )
    ),

    nav_panel(
      "Pathway Enrichment",
      icon = icon("diagram-project"),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("GO Biological Process"),
          card_body(plotOutput("go_plot", height = "500px"))
        ),
        card(
          card_header("KEGG Pathways"),
          card_body(plotOutput("kegg_plot", height = "500px"))
        )
      ),
      card(
        card_header("GO Enrichment Map"),
        card_body(plotOutput("emap_plot", height = "500px"))
      )
    ),

    nav_panel(
      "LLM Interpretation",
      icon = icon("robot"),
      card(
        card_header("Overall Biological Interpretation"),
        card_body(
          div(
            class = "alert alert-warning",
            icon("triangle-exclamation"), " ",
            "AI-generated content. Reviewed for biological accuracy but should be
             critically evaluated in the context of your experimental system."
          ),
          uiOutput("llm_overall")
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Upregulated — Biological Context"),
          card_body(uiOutput("llm_up"))
        ),
        card(
          card_header("Downregulated — Biological Context"),
          card_body(uiOutput("llm_dn"))
        )
      )
    ),

    nav_panel(
      "Log",
      icon = icon("terminal"),
      card(
        card_header("Analysis Log"),
        card_body(verbatimTextOutput("log_output"))
      )
    )
  )
)

# -------------------------------------------------------------------------
# Server
# -------------------------------------------------------------------------

server <- function(input, output, session) {

  # --- Reactive log ---
  log_msgs <- reactiveVal(character(0))

  add_log <- function(msg) {
    ts  <- format(Sys.time(), "[%H:%M:%S]")
    log_msgs(c(log_msgs(), paste(ts, msg)))
  }

  output$log_output <- renderText({
    paste(log_msgs(), collapse = "\n")
  })

  # --- Helper: read CSV/TSV with gzip support ---
  read_input_file <- function(filepath, filename) {
    ext <- tools::file_ext(filename)
    is_gz <- grepl("\\.gz$", filename, ignore.case = TRUE)

    if (is_gz) {
      base_ext <- tools::file_ext(sub("\\.gz$", "", filename, ignore.case = TRUE))
      sep <- if (base_ext == "csv") "," else "\t"
      con <- gzfile(filepath, "rt")
      on.exit(close(con))
      read.table(
        con,
        header = TRUE, sep = sep,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    } else {
      sep <- if (ext == "csv") "," else "\t"
      read.table(
        filepath,
        header = TRUE, sep = sep,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  # --- Read counts ---
  counts_data <- reactive({
    req(input$counts_file)
    read_input_file(input$counts_file$datapath, input$counts_file$name)
  })

  # --- Read metadata ---
  meta_data <- reactive({
    req(input$meta_file)
    read_input_file(input$meta_file$datapath, input$meta_file$name)
  })

  # --- Populate condition selector ---
  output$condition_select <- renderUI({
    req(meta_data())
    cols <- colnames(meta_data())[-1]
    selectInput("condition", "Condition Column", choices = cols)
  })

  # --- Populate reference level selector ---
  output$ref_level_select <- renderUI({
    req(meta_data(), input$condition)
    lvls <- unique(meta_data()[[input$condition]])
    selectInput("ref_level", "Reference Level", choices = lvls)
  })

  # --- Populate model selector ---
  output$model_select <- renderUI({
    models <- list_ollama_models()
    if (length(models) == 0) {
      return(div(
        class = "alert alert-danger p-2",
        icon("circle-xmark"), " Ollama not detected or no models available"
      ))
    }
    selectInput("ollama_model", "Ollama Model", choices = models)
  })

  # --- Main analysis results (reactive) ---
  analysis_results <- reactiveVal(NULL)

  # --- Run analysis on button click ---
  observeEvent(input$run_btn, {

    req(counts_data(), meta_data(), input$condition, input$ref_level)

    analysis_results(NULL)
    log_msgs(character(0))

    withProgress(message = "Running analysis...", value = 0, {

      tryCatch({

        # 1. Validate
        add_log("Validating inputs...")
        incProgress(0.05)

        validation <- validate_inputs(counts_data(), meta_data(), input$condition)
        if (!validation$valid) stop(validation$msg)
        add_log("Inputs validated.")

        # 2. Prepare DESeq2 dataset
        add_log("Preparing DESeq2 dataset...")
        incProgress(0.1)

        dds <- prepare_dds(
          counts_data(), meta_data(),
          input$condition, input$ref_level
        )
        add_log(paste("Genes after pre-filtering:", nrow(dds)))

        # 3. Run DESeq2
        add_log("Running DESeq2...")
        incProgress(0.2)

        deseq_res <- run_deseq2(
          dds, input$condition, input$ref_level,
          input$padj_cutoff, input$lfc_cutoff
        )
        add_log(paste(
          "DESeq2 complete.",
          nrow(deseq_res$significant), "significant genes found."
        ))

        # 4. Enrichment
        add_log("Running GO enrichment...")
        incProgress(0.4)

        all_genes <- deseq_res$results$gene
        sig_genes <- deseq_res$significant$gene

        ego <- run_go_enrichment(sig_genes, all_genes, input$padj_cutoff)
        if (is.null(ego)) add_log("No significant GO terms found.") else
          add_log(paste(nrow(ego@result), "GO terms enriched."))

        add_log("Running KEGG enrichment...")
        incProgress(0.5)

        ekegg <- run_kegg_enrichment(sig_genes, all_genes, input$padj_cutoff)
        if (is.null(ekegg)) add_log("No significant KEGG pathways found.") else
          add_log(paste(nrow(ekegg@result), "KEGG pathways enriched."))

        # 5. Build prompts
        add_log("Building LLM prompts...")
        incProgress(0.6)

        pathway_terms <- extract_top_terms(ego, ekegg)

        prompt_overall <- build_interpretation_prompt(deseq_res, pathway_terms)
        prompt_up      <- build_upregulated_prompt(deseq_res)
        prompt_dn      <- build_downregulated_prompt(deseq_res)

        # 6. Query Ollama (parallel via future/promises)
        add_log(paste("Querying Ollama model:", input$ollama_model, "(3 prompts in parallel)"))
        incProgress(0.65)

        llm_futures <- parallel_query_ollama(
          prompts     = list(prompt_overall, prompt_up, prompt_dn),
          model       = input$ollama_model,
          temperature = input$temperature
        )

        llm_overall <- llm_futures[[1]]
        llm_up      <- llm_futures[[2]]
        llm_dn      <- llm_futures[[3]]

        add_log("All LLM interpretations received.")
        incProgress(0.95)

        # 7. Store results
        analysis_results(list(
          deseq_res   = deseq_res,
          ego         = ego,
          ekegg       = ekegg,
          llm_overall = llm_overall,
          llm_up      = llm_up,
          llm_dn      = llm_dn,
          model_used  = input$ollama_model
        ))

        add_log("Analysis complete.")
        incProgress(1)

      }, error = function(e) {
        add_log(paste("ERROR:", conditionMessage(e)))
        showNotification(conditionMessage(e), type = "error", duration = 10)
      })
    })
  })

  # --- Plots ---
  output$pca_plot <- renderPlot({
    req(analysis_results())
    plot_pca(analysis_results()$deseq_res$dds, input$condition)
  })

  output$volcano_plot <- renderPlot({
    req(analysis_results())
    r <- analysis_results()$deseq_res
    plot_volcano(
      r$results,
      padj_cutoff = input$volcano_padj,
      lfc_up      = input$volcano_lfc_up,
      lfc_dn      = input$volcano_lfc_dn
    )
  })

  output$go_plot <- renderPlot({
    req(analysis_results())
    plot_go_dotplot(analysis_results()$ego)
  })

  output$kegg_plot <- renderPlot({
    req(analysis_results())
    plot_kegg_dotplot(analysis_results()$ekegg)
  })

  output$emap_plot <- renderPlot({
    req(analysis_results())
    plot_go_emap(analysis_results()$ego)
  })

  # --- Tables (react to volcano cutoff inputs) ---
  output$up_table <- renderTable({
    req(analysis_results())
    analysis_results()$deseq_res$results |>
      dplyr::filter(padj < input$volcano_padj, log2FoldChange > input$volcano_lfc_up) |>
      dplyr::arrange(padj) |>
      dplyr::slice_head(n = 20) |>
      dplyr::mutate(
        log2FoldChange = round(log2FoldChange, 3),
        padj           = signif(padj, 3)
      ) |>
      dplyr::select(gene, log2FoldChange, padj)
  }, striped = TRUE, hover = TRUE)

  output$dn_table <- renderTable({
    req(analysis_results())
    analysis_results()$deseq_res$results |>
      dplyr::filter(padj < input$volcano_padj, log2FoldChange < input$volcano_lfc_dn) |>
      dplyr::arrange(padj) |>
      dplyr::slice_head(n = 20) |>
      dplyr::mutate(
        log2FoldChange = round(log2FoldChange, 3),
        padj           = signif(padj, 3)
      ) |>
      dplyr::select(gene, log2FoldChange, padj)
  }, striped = TRUE, hover = TRUE)

  # --- LLM outputs ---
  output$llm_overall <- renderUI({
    req(analysis_results())
    div(
      style = "white-space: pre-wrap; font-size: 0.95rem; line-height: 1.6;",
      analysis_results()$llm_overall
    )
  })

  output$llm_up <- renderUI({
    req(analysis_results())
    div(
      style = "white-space: pre-wrap; font-size: 0.95rem; line-height: 1.6;",
      analysis_results()$llm_up
    )
  })

  output$llm_dn <- renderUI({
    req(analysis_results())
    div(
      style = "white-space: pre-wrap; font-size: 0.95rem; line-height: 1.6;",
      analysis_results()$llm_dn
    )
  })

  # --- Download button (visible only after analysis completes) ---
  output$download_ui <- renderUI({
    req(analysis_results())
    downloadButton("dl_html", "Download HTML Report", class = "btn-success w-100")
  })

  output$dl_html <- downloadHandler(
    filename = function() paste0("rnaseq_report_", Sys.Date(), ".html"),
    content  = function(file) {
      req(analysis_results())
      r <- analysis_results()

      report_params <- list(
        results            = r$deseq_res$results,
        significant        = r$deseq_res$significant,
        contrast_level     = r$deseq_res$contrast_level,
        ref_level          = r$deseq_res$ref_level,
        padj_cutoff        = r$deseq_res$padj_cutoff,
        lfc_cutoff         = r$deseq_res$lfc_cutoff,
        go_results         = if (!is.null(r$ego))   r$ego@result   else NULL,
        kegg_results       = if (!is.null(r$ekegg)) r$ekegg@result else NULL,
        llm_interpretation = r$llm_overall,
        llm_upregulated    = r$llm_up,
        llm_downregulated  = r$llm_dn,
        model_used         = r$model_used
      )

      tmp_dir <- tempdir()
      tmp_rmd <- file.path(tmp_dir, "report.Rmd")
      file.copy(here("report", "report.Rmd"), tmp_rmd, overwrite = TRUE)

      rmarkdown::render(
        input         = tmp_rmd,
        output_format = "html_document",
        output_file   = file.path(tmp_dir, "report.html"),
        params        = report_params,
        envir         = new.env(parent = globalenv()),
        quiet         = TRUE
      )

      file.copy(file.path(tmp_dir, "report.html"), file)
    }
  )

}

# -------------------------------------------------------------------------
shinyApp(ui, server)
