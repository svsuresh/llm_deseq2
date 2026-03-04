# ollama.R
# Bridge between R and local Ollama instance via httr2

library(httr2)

#' Check if Ollama is running and model is available
#' @param model Character. Model name to check
#' @return Logical
check_ollama <- function(model = "mistral") {
  tryCatch({
    resp <- request("http://localhost:11434/api/tags") |>
      req_perform() |>
      resp_body_json()
    
    available <- sapply(resp$models, function(m) m$name)
    any(grepl(model, available))
  }, error = function(e) {
    FALSE
  })
}

#' Query local Ollama model
#' @param prompt Character. The prompt to send
#' @param model Character. Ollama model name
#' @param temperature Numeric. Sampling temperature (0-1)
#' @return Character. Model response text
query_ollama <- function(prompt, 
                         model = "mistral",
                         temperature = 0.3) {
  tryCatch({
    resp <- request("http://localhost:11434/api/chat") |>
      req_body_json(list(
        model = model,
        messages = list(
          list(
            role = "system",
            content = paste(
              "You are an expert bioinformatics analyst.",
              "Provide concise, accurate biological interpretations.",
              "Do not speculate beyond the data provided.",
              "Do not fabricate gene functions or pathway associations.",
              "Use precise scientific language appropriate for a research audience."
            )
          ),
          list(
            role = "user",
            content = prompt
          )
        ),
        options = list(temperature = temperature),
        stream = FALSE
      )) |>
      req_timeout(120) |>
      req_perform() |>
      resp_body_json()
    
    resp$message$content
    
  }, error = function(e) {
    stop(paste("Ollama query failed:", conditionMessage(e)))
  })
}

#' Get list of available models from Ollama
#' @return Character vector of model names
list_ollama_models <- function() {
  tryCatch({
    resp <- request("http://localhost:11434/api/tags") |>
      req_perform() |>
      resp_body_json()
    
    sapply(resp$models, function(m) m$name)
  }, error = function(e) {
    character(0)
  })
}
