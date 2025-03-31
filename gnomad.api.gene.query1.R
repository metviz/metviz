# Load required libraries
if (!require("httr")) install.packages("httr")
if (!require("jsonlite")) install.packages("jsonlite")
if (!require("stringr")) install.packages("stringr")

library(httr)
library(jsonlite)
library(stringr)

# Read the GraphQL query template from a file
read_query_template <- function(template_file) {
  if (!file.exists(template_file)) {
    stop("Template file does not exist")
  }
  query <- readLines(template_file, warn = FALSE)
  return(paste(query, collapse = "\n"))
}

# Replace the GeneSymbol placeholder in the query
prepare_query <- function(query_template, gene_symbol) {
  if (missing(gene_symbol) || gene_symbol == "") {
    stop("GeneSymbol must be provided")
  }
  query <- str_replace_all(query_template, fixed("{{GeneSymbol}}"), gene_symbol)
  return(query)
}

# Main function to execute the API query
query_gnomad_api <- function(query) {
  url <- "https://gnomad.broadinstitute.org/api"
  response <- POST(
    url,
    content_type("application/graphql; charset=utf-8"),
    body = query,
    encode = "raw"
  )
  return(response)
}

# Main execution
tryCatch({
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 2) {
    cat("Enter the path to the GraphQL query template file: ")
    template_file <- readLines(con = stdin(), n = 1)
    
    cat("Enter the GeneSymbol (e.g., BRCA1): ")
    gene_symbol <- readLines(con = stdin(), n = 1)
  } else {
    template_file <- args[1]
    gene_symbol <- args[2]
  }

  json_file <- paste0(gene_symbol, ".gnomad_variants_data.json")

  if (file.exists(json_file)) {
    cat("Using existing JSON file...\n")
    result <- fromJSON(json_file)
  } else {
    query_template <- read_query_template(template_file)
    query <- prepare_query(query_template, gene_symbol)

    response <- query_gnomad_api(query)
    if (status_code(response) != 200) {
      stop(paste("API Error:", content(response, "text")))
    }

    result <- fromJSON(content(response, "text"))
    write(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE), json_file)
  }

  cat("Raw JSON response saved to:", json_file, "\n")
}, error = function(e) {
  cat("Error occurred:", conditionMessage(e), "\n")
})

