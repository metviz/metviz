library(ggplot2)
library(plotly)
library(readr)
library(dplyr)

# Function to get mean pathogenicity data
get_mean_pathogenicity <- function(uniprot_id) {
  url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", uniprot_id, "-F1-aa-substitutions.csv")
  afs_data <- read_csv(url)
  
  # Compute mean pathogenicity per position
  mean_data <- afs_data %>%
    mutate(Position = as.numeric(substr(protein_variant, 2, nchar(protein_variant)-1))) %>%
    group_by(Position) %>%
    summarise(mean_pathogenicity = mean(am_pathogenicity, na.rm = TRUE),
              variants = paste(protein_variant[am_pathogenicity >= 0.564], collapse = ", ")) %>%
    mutate(
      category = case_when(
        mean_pathogenicity < 0.34 ~ "Benign",
        mean_pathogenicity >= 0.34 & mean_pathogenicity < 0.564 ~ "Uncertain",
        mean_pathogenicity >= 0.564 ~ "Pathogenic"
      )
    )
  
  return(mean_data)
}

# Example UniProt ID
uniprot_id <- "Q8NBP7"
mean_data <- get_mean_pathogenicity(uniprot_id)

# Define color scale matching the image
color_scale <- c(
  "Benign" = "#1E90FF",      # Dodger blue
  "Uncertain" = "#FFD700",   # Gold
  "Pathogenic" = "#FF4500"   # OrangeRed
)

# Create bar chart with gradient coloring
afmpsgg <- ggplot(mean_data, aes(x = Position, y = mean_pathogenicity, fill = category,
                            text = paste("Position:", Position,
                                         "<br>Mean Pathogenicity:", round(mean_pathogenicity, 3),
                                         "<br>Category:", category,
                                         "<br>Pathogenic Variants:", variants))) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = color_scale) +
  geom_hline(yintercept = 0.34, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0.564, linetype = "dashed", color = "gray50") +
  #annotate("text", x = max(mean_data$Position)*0.1, y = 0.17, label = "Benign", color = "gray30") +
  #annotate("text", x = max(mean_data$Position)*0.1, y = 0.45, label = "Uncertain", color = "gray30") +
  #annotate("text", x = max(mean_data$Position)*0.1, y = 0.78, label = "Pathogenic", color = "gray30") +
  labs(x = "Residue Position", y = "Mean Pathogenicity",
       title = paste("Mean Pathogenicity Scores for", uniprot_id)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

# Convert to interactive plot
plotly_afmps <- ggplotly(afmpsgg, tooltip = "text") %>%
  layout(
    hoverlabel = list(bgcolor = "white"),
    xaxis = list(tickmode = "linear"),
    shapes = list(
      list(type = "line", y0 = 0.34, y1 = 0.34, x0 = 0, x1 = 1, xref = "paper",
           line = list(color = "gray", dash = "dash")),
      list(type = "line", y0 = 0.564, y1 = 0.564, x0 = 0, x1 = 1, xref = "paper",
           line = list(color = "gray", dash = "dash"))
    ),
    annotations = list(
      list(x = 0.05, y = 0.17, xref = "paper", yref = "y",
           text = "Benign", showarrow = FALSE, font = list(color = "gray30")),
      list(x = 0.05, y = 0.45, xref = "paper", yref = "y",
           text = "Uncertain", showarrow = FALSE, font = list(color = "gray30")),
      list(x = 0.05, y = 0.78, xref = "paper", yref = "y",
           text = "Pathogenic", showarrow = FALSE, font = list(color = "gray30"))
    )
  )

# Display plot
plotly_afmps
