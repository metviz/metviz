library(httr)
library(jsonlite)

# Function to fetch UniProt features data for an input UniProt ID
# Available fields from Uniprot
# https://www.uniprot.org/help/return_fields
# https://www.uniprot.org/uniprotkb/Q8NBP7/entry

extract_pfam <- function(x) { 
  tryCatch({
    resp <- request("https://rest.uniprot.org/uniprotkb/") %>%
      req_url_path_append(paste0(x, ".json")) %>%
      req_url_query(fields = "accession,id,protein_name,gene_primary,ft_repeat,ft_region,ft_domain,ft_compbias,ft_coiled,ft_motif,ft_zn_fing,length") %>%
      req_error(body = function(resp) {
        paste("Error:", resp_status(resp), resp_status_desc(resp))
      }) %>%
      req_perform()
    
    # Parse Uniprot JSON response
    pfam_data <- resp %>%
      resp_body_string() %>%
      fromJSON()
    
    # Define pfam domain color palette 
    color_palette <- c("#2dcf00", "#ff5353", "#5b5bff", "#ebd61d", "#ba21e0",
                       "#ff9c42", "#ff7dff", "#b9264f", "#baba21", "#c48484",
                       "#1f88a7", "#cafeb8", "#4a9586", "#ceb86c", "#0e180f")
    # Create a list to store used colors for each feature type
    used_colors <- list()
    
    # Function to process uniprot features types 
    process_feature_type <- function(features, feature_type) {
      type_features <- features[features$type == feature_type, ]
      if (nrow(type_features) == 0) return(NULL)
      
      # Get available colors (excluding those used by other types)
      available_colors <- setdiff(color_palette, unlist(used_colors))
      
      # If we've used all colors, reset the available colors
      if (length(available_colors) == 0) {
        available_colors <- color_palette
      }
      
      # Assign colors, cycling through the available colors if needed
      colors <- rep_len(available_colors, nrow(type_features))
      
      # Used colors for this feature type
      used_colors[[feature_type]] <<- unique(colors)
      
      data.frame(
        type = type_features$type,
        description = type_features$description,
        start = type_features$location$start$value,
        end = type_features$location$end$value,
        color = colors,
        stringsAsFactors = FALSE
      )
    }
    
    # Get unique feature types
    feature_types <- unique(pfam_data$features$type)
    
    # Process each feature type and add to pfam_data
    for (feature_type in feature_types) {
      processed_df <- process_feature_type(pfam_data$features, feature_type)
      if (!is.null(processed_df)) {
        pfam_data[[feature_type]] <- processed_df
      }
    }
    
    # Protein Length
    pfam_data$length <- pfam_data$sequence$length
    return(pfam_data)
    
  }, error = function(e) {
    stop(paste("Request failed:", e$message))
  })
}


# Extract protein features data for UniProt ID P43694
pfam_data <- extract_pfam("Q8NBP7")

# Access the processed data frames:
print(pfam_data$uniProtkbId)
print(pfam_data$primaryAccession)
print(pfam_data$proteinDescription$recommendedName$fullName$value)
print(pfam_data$length)

# Feature info
print(pfam_data$Domain)
print(pfam_data$Region)
print(pfam_data$`Zinc finger`)
print(pfam_data$`Compositional bias`)
