library(jsonlite)
library(dplyr)
library(stringr)
library(purrr)

# Update Working Directory
setwd("/home/metviz/gnomad")  

# 
parse_gnomad_json <- function(json_file) {
  # 1. Read and clean JSON file
  json_text <- tryCatch({
    readLines(json_file, warn = FALSE) %>% 
      paste(collapse = "\n") %>%
      str_replace_all("'", "\"") %>%
      str_replace_all("None", "null") %>%
      str_replace_all("True", "true") %>%
      str_replace_all("False", "false") %>%
      str_replace_all(",\\s*\\}", "}") %>%
      str_replace_all(",\\s*\\]", "]") %>%
      str_replace_all(":\\s*NaN", ": null")
  }, error = function(e) {
    stop(paste("Failed to read JSON file:", e$message))
  })
  
  # 2. Parse JSON with structure validation
  data <- tryCatch({
    parsed <- fromJSON(json_text, simplifyDataFrame = FALSE)
    if (is.null(parsed$data$gene$variants)) {
      stop("No variants found in the JSON structure")
    }
    parsed$data$gene$variants
  }, error = function(e) {
    stop(paste("JSON parsing failed:", e$message))
  })
  
  # 3. Safe extraction function
  safe_extract <- function(x, field, default = NA, type = NULL) {
    tryCatch({
      val <- x[[field]]
      if (is.null(val)) return(default)
      if (identical(val, list())) return(default)
      
      if (!is.null(type)) {
        if (type == "int") val <- as.integer(val)
        if (type == "num") val <- as.numeric(val)
        if (type == "char") val <- as.character(val)
        if (type == "flag") val <- paste(unlist(val), collapse = ";")
      }
      val
    }, error = function(e) default)
  }
  
  # 4. Process each variant with all requested fields
  variants_df <- map_df(data, function(variant) {
    # Extract basic variant info
    variant_parts <- str_split_fixed(safe_extract(variant, "variant_id", "", "char"), "-", 4)
    
    # Process population data
    process_populations <- function(pop_data) {
      if (!is.list(pop_data)) return(NA_character_)
      pops <- map_chr(pop_data, function(pop) {
        paste(safe_extract(pop, "id", "unknown", "char"),
              safe_extract(pop, "ac_hom", 0, "int"),
              sep = ":")
      })
      paste(pops, collapse = ",")
    }
    
    # Process in-silico predictors
    process_predictors <- function(pred_data) {
      if (!is.list(pred_data)) return(list())
      setNames(
        map_chr(pred_data, ~ safe_extract(.x, "value", NA, "char")),
        map_chr(pred_data, ~ safe_extract(.x, "id", "unknown", "char"))
      )
    }
    
    predictors <- process_predictors(safe_extract(variant, "in_silico_predictors", list()))
    
    # Create data frame with all requested fields
    data.frame(
      # Basic variant info
      variant_id = safe_extract(variant, "variant_id", NA, "char"),
      chrom = "1",
      position = safe_extract(variant, "pos", NA, "int"),
      ref = ifelse(ncol(variant_parts) >= 3, variant_parts[3], NA),
      alt = ifelse(ncol(variant_parts) >= 4, variant_parts[4], NA),
      
      # Transcript info
      transcript_id = safe_extract(variant, "transcript_id", NA, "char"),
      hgvsc = safe_extract(variant, "hgvsc", NA, "char"),
      hgvsp = safe_extract(variant, "hgvsp", NA, "char"),
      consequence = safe_extract(variant, "consequence", NA, "char"),
      flags = safe_extract(variant, "flags", NA, "flag"),
      
      # Exome data
      exome_ac = safe_extract(variant$exome, "ac", 0, "int"),
      exome_an = safe_extract(variant$exome, "an", 0, "int"),
      exome_af = safe_extract(variant$exome, "af", 0, "num"),
      exome_ac_hom = safe_extract(variant$exome, "ac_hom", 0, "int"),
      exome_filters = paste(safe_extract(variant$exome, "filters", list(), "list"), collapse = ";"),
      
      # Population data
      populations_ac_hom = process_populations(safe_extract(variant$exome, "populations", list())),
      
      # In-silico predictors
      cadd_score = safe_extract(predictors, "cadd", NA, "num"),
      spliceai_score = safe_extract(predictors, "spliceai_ds_max", NA, "num"),
      pangolin_score = safe_extract(predictors, "pangolin_largest_ds", NA, "num"),
      phylop_score = safe_extract(predictors, "phylop", NA, "num"),
      
      stringsAsFactors = FALSE
    )
  })
  
  # 5. Post-processing
  # Replace NA with 0 for count columns
  count_cols <- c("exome_ac", "exome_an", "exome_ac_hom")
  for (col in count_cols) {
    variants_df[[col]][is.na(variants_df[[col]])] <- 0
  }
  
  # Replace NA with 0 for frequency columns
  freq_cols <- c("exome_af")
  for (col in freq_cols) {
    variants_df[[col]][is.na(variants_df[[col]])] <- 0
  }
  
  # Remove duplicate variants
  variants_df <- distinct(variants_df, variant_id, .keep_all = TRUE)
  
  # Add metadata
  attr(variants_df, "extracted_on") <- Sys.time()
  attr(variants_df, "source_file") <- json_file
  attr(variants_df, "total_variants") <- nrow(variants_df)
  
  return(variants_df)
}

# Run the parser with error handling
tryCatch({
  cat("Starting gnomAD JSON parsing...\n")
  start_time <- Sys.time()
  
  ## INPUT Gene GnomAD JSON file (downloaded using API) ##
  variants_data <- parse_gnomad_json("PCSK9.gnomad_variants_data.json")
  
  cat(sprintf("\nSuccessfully processed %d variants (%.1f seconds)\n", 
              nrow(variants_data), 
              difftime(Sys.time(), start_time, units = "secs")))
  
  cat("\nFirst 3 variants:\n")
  print(head(variants_data, 3))
  
  ## OUTPUT varaints table in comma separated (csv) file ##
  output_file <- paste0("PCSK9_gnomad_variants_", format(Sys.time(), "%Y%m%d"), ".csv")
  write.csv(variants_data, output_file, row.names = FALSE)
  cat(sprintf("\nResults saved to: %s\n", output_file))
  
}, error = function(e) {
  cat(paste("\nERROR:", e$message, "\n"))
  if (exists("variants_data")) {
    cat("\nPartial data structure:\n")
    str(variants_data, max.level = 1)
  }
})
