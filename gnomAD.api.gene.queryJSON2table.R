suppressMessages(library(jsonlite))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(purrr))

parse_gnomad_json <- function(json_file) {
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
  
  data <- tryCatch({
    parsed <- fromJSON(json_text, simplifyDataFrame = FALSE)
    if (is.null(parsed$data$gene$variants)) {
      stop("No variants found in the JSON structure")
    }
    parsed$data$gene$variants
  }, error = function(e) {
    stop(paste("JSON parsing failed:", e$message))
  })
  
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
  
  variants_df <- map_df(data, function(variant) {
    variant_parts <- str_split_fixed(safe_extract(variant, "variant_id", "", "char"), "-", 4)
    
    process_populations <- function(pop_data) {
      if (!is.list(pop_data)) return(NA_character_)
      pops <- map_chr(pop_data, function(pop) {
        paste(safe_extract(pop, "id", "unknown", "char"),
              safe_extract(pop, "ac_hom", 0, "int"),
              sep = ":")
      })
      paste(pops, collapse = ",")
    }
    
    process_predictors <- function(pred_data) {
      if (!is.list(pred_data)) return(list())
      setNames(
        map_chr(pred_data, ~ safe_extract(.x, "value", NA, "char")),
        map_chr(pred_data, ~ safe_extract(.x, "id", "unknown", "char"))
      )
    }
    
    predictors <- process_predictors(safe_extract(variant, "in_silico_predictors", list()))
    
    data.frame(
      variant_id = safe_extract(variant, "variant_id", NA, "char"),
      chrom = "1",
      position = safe_extract(variant, "pos", NA, "int"),
      ref = ifelse(ncol(variant_parts) >= 3, variant_parts[3], NA),
      alt = ifelse(ncol(variant_parts) >= 4, variant_parts[4], NA),
      
      transcript_id = safe_extract(variant, "transcript_id", NA, "char"),
      hgvsc = safe_extract(variant, "hgvsc", NA, "char"),
      hgvsp = safe_extract(variant, "hgvsp", NA, "char"),
      consequence = safe_extract(variant, "consequence", NA, "char"),
      flags = safe_extract(variant, "flags", NA, "flag"),
      
      exome_ac = safe_extract(variant$exome, "ac", 0, "int"),
      exome_an = safe_extract(variant$exome, "an", 0, "int"),
      exome_af = safe_extract(variant$exome, "af", 0, "num"),
      exome_ac_hom = safe_extract(variant$exome, "ac_hom", 0, "int"),
      exome_filters = paste(safe_extract(variant$exome, "filters", list(), "list"), collapse = ";"),
      
      populations_ac_hom = process_populations(safe_extract(variant$exome, "populations", list())),
      
      cadd_score = safe_extract(predictors, "cadd", NA, "num"),
      spliceai_score = safe_extract(predictors, "spliceai_ds_max", NA, "num"),
      pangolin_score = safe_extract(predictors, "pangolin_largest_ds", NA, "num"),
      phylop_score = safe_extract(predictors, "phylop", NA, "num"),
      
      stringsAsFactors = FALSE
    )
  })
  
  # Filter for hgvsp entries
  variants_df <- variants_df %>% 
    filter(!is.na(hgvsp)) %>%
    distinct(variant_id, .keep_all = TRUE)
  
  # Parse hgvsp for aa_ref, aa_pos, aa_alt
  three_to_one <- c(
    Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
    Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
    Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
    Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
    Ter = "*", Sec = "U"
  )
  
  parse_hgvsp <- function(hgvsp_string) {
    hgvsp_string <- str_extract(hgvsp_string, "p\\.[^ ]+")
    
    if (is.na(hgvsp_string)) return(c(NA, NA, NA))
    
    # 3-letter to 1-letter amino acid code
    three_to_one <- c(
      Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
      Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
      Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
      Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
      Ter = "*", Sec = "U"
    )
    
    # Match: p.3letter###3letter or p.3letter###X or other complex types
    if (str_detect(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)([A-Za-z]{3}|X|Ter).*")) {
      parts <- str_match(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)([A-Za-z]{3}|X|Ter)")
      aa_ref <- three_to_one[parts[2]]
      aa_pos <- as.integer(parts[3])
      aa_alt <- ifelse(parts[4] == "X" || parts[4] == "Ter", "*", three_to_one[parts[4]])
      return(c(aa_ref, aa_pos, aa_alt))
    }
    
    # Deletions like p.Ile507del or p.Phe508del
    if (str_detect(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)del")) {
      parts <- str_match(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)del")
      aa_ref <- three_to_one[parts[2]]
      aa_pos <- as.integer(parts[3])
      return(c(aa_ref, aa_pos, "del"))
    }
    
    # Delins like p.Ser607delinsTer
    if (str_detect(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)delins([A-Za-z]{3}|Ter|X)")) {
      parts <- str_match(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)delins([A-Za-z]{3}|Ter|X)")
      aa_ref <- three_to_one[parts[2]]
      aa_pos <- as.integer(parts[3])
      aa_alt <- ifelse(parts[4] %in% c("Ter", "X"), "*", three_to_one[parts[4]])
      return(c(aa_ref, aa_pos, aa_alt))
    }
    
    # Frameshifts like p.Glu1011AlafsTer14 or p.Lys684AsnfsX38
    if (str_detect(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)[A-Za-z]{3}?fs")) {
      parts <- str_match(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)[A-Za-z]{3}?fs")
      aa_ref <- three_to_one[parts[2]]
      aa_pos <- as.integer(parts[3])
      return(c(aa_ref, aa_pos, "fs"))
    }
    
    # Start loss: p.Met1?
    if (str_detect(hgvsp_string, "p\\.Met1\\?")) {
      return(c("M", 1, "?"))
    }
    
    # Unknown or unsupported format
    return(c(NA, NA, NA))
  }
  
  
  parsed_aas <- t(sapply(variants_df$hgvsp, parse_hgvsp))
  variants_df$aa_ref <- parsed_aas[, 1]
  variants_df$aa_pos <- as.integer(parsed_aas[, 2])
  variants_df$aa_alt <- parsed_aas[, 3]
  
  # Replace NA with 0 for counts and frequencies
  count_cols <- c("exome_ac", "exome_an", "exome_ac_hom")
  freq_cols <- c("exome_af")
  for (col in count_cols) {
    variants_df[[col]][is.na(variants_df[[col]])] <- 0
  }
  for (col in freq_cols) {
    variants_df[[col]][is.na(variants_df[[col]])] <- 0
  }
  
  attr(variants_df, "extracted_on") <- Sys.time()
  attr(variants_df, "source_file") <- json_file
  attr(variants_df, "total_variants") <- nrow(variants_df)
  
  return(variants_df)
}

# --------- MAIN EXECUTION ---------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript gnomad.api.gene.query2table.R <Gene.gnomad_variants_data.json>")
}

input_file <- args[1]
gene_symbol <- str_replace(basename(input_file), "\\.gnomad_variants_data\\.json$", "")
output_file <- paste0(gene_symbol, "_gnomad_variants_with_hgvsp_", format(Sys.Date(), "%Y%m%d"), ".csv")

tryCatch({
  cat("Parsing file:", input_file, "\n")
  start_time <- Sys.time()
  
  variants_data <- parse_gnomad_json(input_file)
  
  cat(sprintf("Processed %d variants with hgvsp (%.1f seconds)\n", 
              nrow(variants_data), 
              difftime(Sys.time(), start_time, units = "secs")))
  
  write.csv(variants_data, output_file, row.names = FALSE)
  cat("Saved CSV to:", output_file, "\n")
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
})
