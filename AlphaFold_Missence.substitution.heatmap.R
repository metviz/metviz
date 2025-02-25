library(ggplot2)
library(plotly)
library(readr)

# Function to read data from URL based on UniProt ID
get_protein_data <- function(uniprot_id) {
  url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", uniprot_id, "-F1-aa-substitutions.csv")
  data <- read_csv(url)
  return(data)
}

# Example UniProt ID
uniprot_id <- "Q9NZC2"
data <- get_protein_data(uniprot_id)

# Split 'protein_variant' column into Ref AA, Position, and Var AA
data$Ref_AA <- substr(data$protein_variant, 1, 1)
data$Position <- as.numeric(substr(data$protein_variant, 2, nchar(data$protein_variant) - 1))
data$Var_AA <- substr(data$protein_variant, nchar(data$protein_variant), nchar(data$protein_variant))

# Create heatmap plot
gg <- ggplot(data, aes(x = factor(Position), y = Var_AA, fill = am_pathogenicity, text = paste("Variant:", protein_variant, "\nPathogenicity:", am_pathogenicity, "\nClass:", am_class))) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Residue Position", y = "Variant Amino Acid", fill = "Pathogenicity", title = "Pathogenicity Heatmap per Position") +
  theme_minimal()

# Convert to interactive plotly plot
plotly_gg <- ggplotly(gg, tooltip = "text")

# Display plot
plotly_gg
