library(tidyverse)
library(plotly)

# Function to read data from AlphaFold Protein Structure Database hosted at EBI and calculate mean pathogenicity per residue
extract_af_missense_path <- function(uniprot_id) {
  url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", uniprot_id, "-F1-aa-substitutions.csv")
  
  # Read data safely
  afaasub <- tryCatch({
    read_csv(url, show_col_types = FALSE)
  }, error = function(e) {
    message("Error fetching data: ", e$message)
    return(NULL)
  })
  
  if (is.null(afaasub)) return(NULL)
  
  # Extract amino acid change details
  afaasub <- afaasub %>%
    mutate(
      Ref_AA = substr(protein_variant, 1, 1),
      Position = as.numeric(substr(protein_variant, 2, nchar(protein_variant) - 1)),
      Var_AA = substr(protein_variant, nchar(protein_variant), nchar(protein_variant))
    )
  
  # Compute mean pathogenicity per position
  mean_pathogenicity <- afaasub %>%
    group_by(Position) %>%
    summarise(
      am_pathogenicity = mean(am_pathogenicity, na.rm = TRUE),  # Mean pathogenicity
      protein_variant = paste(protein_variant[am_pathogenicity >= 0.564], collapse = ", ")  # Keep high pathogenic variants
    ) %>%
    mutate(Var_AA = "Mean")  # Rename 1Mean to Mean
  
  return(mean_pathogenicity)  # Return only mean pathogenicity data
}

# Example usage
uniprot_id <- "Q9NZC2"
afaasub_mean <- extract_af_missense_path(uniprot_id)

# Create heatmap plot (Mean row only)
gg_mean <- ggplot(afaasub_mean, aes(
  x = factor(Position), y = Var_AA, fill = am_pathogenicity, 
  text = paste0("<b>Position:</b> ", Position,
                "<br><b>Mean Pathogenicity:</b> ", round(am_pathogenicity, 3), 
                "<br><b>AlphaMissense Pathogenic Variants:</b> ", protein_variant)
)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "gray") +
  labs(
    x = "Residue Position", 
    y = "Mean Pathogenicity", 
    fill = "Pathogenicity", 
    title = "Mean Pathogenicity Heatmap"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.ticks.x = element_blank())

# Convert to interactive plotly plot
plotly_gg_mean <- ggplotly(gg_mean, tooltip = "text")

# Display plot
plotly_gg_mean
