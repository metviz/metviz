library(ggplot2)
library(plotly)
library(readr)
library(dplyr)

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

# Compute mean pathogenicity per position
mean_pathogenicity <- data %>%
  group_by(Position) %>%
  summarise(am_pathogenicity = mean(am_pathogenicity, na.rm = TRUE),
            protein_variant = paste(protein_variant[am_pathogenicity >= 0.564], collapse=", ")) %>%
  mutate(Var_AA = "1Mean") # mean represented as 1mean to display mean as the bottom row

# Combine mean pathogenicity row with original data
data <- bind_rows(mean_pathogenicity, data) %>% arrange(Position, desc(Var_AA))

# Create heatmap plot
gg <- ggplot(data, aes(x = factor(Position), y = Var_AA, fill = am_pathogenicity, 
                       text = ifelse(Var_AA == "1Mean", paste("Mean Pathogenicity:", am_pathogenicity, "\nVariants:", protein_variant), 
                                     paste("Variant:", protein_variant, "\nPathogenicity:", am_pathogenicity, "\nClass:", am_class)))) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Residue Position", y = "Variant Amino Acid", fill = "Pathogenicity", title = "Pathogenicity Heatmap per Position") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Convert to interactive plotly plot
plotly_gg <- ggplotly(gg, tooltip = "text")

# Display plot
plotly_gg
