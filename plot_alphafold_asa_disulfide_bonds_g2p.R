## R script to plot AlphaFold pLDDT, Accessible Surface Area, and Disulfide bonds of a protein data from G2P 
# Load necessary libraries
library(ggplot2)
library(plotly)
library(dplyr)

# Function to extract G2P data
extract_g2p <- function(g, u) {
  url <- paste0("https://g2p.broadinstitute.org/api/gene/", g, "/protein/", u, "/protein-features")
  print(url)
  
  g2p_data <- read.table(url, sep = '\t', header = TRUE, fileEncoding = "utf-8")
  return(g2p_data)
}

# Fetch data
g2p_data <- extract_g2p("PCSK9", "Q8NBP7")

# Extract specific columns for each feature
features <- list(
  Disulfide_bond = if ("Disulfide.bond..UniProt." %in% colnames(g2p_data)) {
    g2p_data[g2p_data$Disulfide.bond..UniProt. != "", c("residueId", "AA", "Disulfide.bond..UniProt.")]
  } else {
    data.frame()
  }
)

# Parse disulfide bond pairs if the column exists
if (nrow(features$Disulfide_bond) > 0) {
  disulfide_pairs <- features$Disulfide_bond %>%
    group_by(Disulfide.bond..UniProt.) %>%
    summarise(residueId1 = min(residueId),
              residueId2 = max(residueId)) %>%
    ungroup() %>%
    distinct(residueId1, residueId2, .keep_all = TRUE)
} else {
  disulfide_pairs <- data.frame()
}

# Define function to classify pLDDT confidence
classify_confidence <- function(pLDDT) {
  case_when(
    pLDDT >= 90 ~ "Very high (90-100)",
    pLDDT >= 70 ~ "Confident (70-90)",
    pLDDT >= 50 ~ "Low (50-70)",
    TRUE ~ "Very low (0-50)"
  )
}

# Correct color scheme based on user specification
confidence_colors <- c(
  "Very low (0-50)" = "#FF6B6B",  # Red
  "Low (50-70)" = "#FFD93D",      # Yellow
  "Confident (70-90)" = "#4D96FF", # Blue
  "Very high (90-100)" = "#6BCB77" # Green
)

# Add classification column with ordered factors
g2p_data$Confidence_Level <- factor(
  classify_confidence(g2p_data$AlphaFold.confidence..pLDDT.),
  levels = names(confidence_colors)
)

# 📌 Function to plot AlphaFold pLDDT
plot_pLDDT <- function(g2p_data) {
  p <- ggplot(g2p_data, aes(x = residueId, y = AlphaFold.confidence..pLDDT., 
                            fill = Confidence_Level,
                            text = paste("Residue:", residueId,
                                         "<br>AA:", AA,
                                         "<br>pLDDT:", round(AlphaFold.confidence..pLDDT., 1),
                                         "<br>Confidence:", Confidence_Level))) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = confidence_colors) +
    labs(title = "AlphaFold Confidence (pLDDT)",
         x = "Residue Position",
         y = "pLDDT Score",
         fill = "Confidence Level") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    scale_x_continuous(expand = expansion(mult = 0.01))
  
  return(ggplotly(p, tooltip = "text") %>%
           layout(hoverlabel = list(bgcolor = "white"),
                  legend = list(orientation = "h", y = -0.2)))
}

# 📌 Function to plot ASA
plot_ASA <- function(g2p_data) {
  p <- ggplot(g2p_data, aes(x = residueId, y = Accessible.surface.area..Å..., 
                            text = paste("Residue:", residueId, 
                                         "<br>AA:", AA,
                                         "<br>ASA:", round(Accessible.surface.area..Å..., 1)))) +
    geom_bar(stat = "identity", fill = "#2ca02c", width = 1) +
    labs(title = "Accessible Surface Area (ASA) Distribution",
         x = "Residue Position",
         y = "ASA (Å²)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12)
    ) +
    scale_y_continuous(expand = expansion(mult = 0.02))
  
  return(ggplotly(p, tooltip = "text"))
}

# 📌 Function to plot Disulfide bonds
plot_disulfide_bonds <- function(disulfide_pairs) {
  if (nrow(disulfide_pairs) == 0) {
    return(NULL)
  }
  
  # Assign a unique y-level for each bond
  disulfide_pairs$y_level <- seq(nrow(disulfide_pairs), 1, by = -1)
  
  p <- ggplot() +
    geom_segment(data = disulfide_pairs, 
                 aes(x = residueId1, y = y_level, xend = residueId2, yend = y_level), 
                 color = "blue", size = 1.2) +
    geom_point(data = disulfide_pairs, 
               aes(x = residueId1, y = y_level), 
               color = "red", size = 3) +
    geom_point(data = disulfide_pairs, 
               aes(x = residueId2, y = y_level), 
               color = "red", size = 3) +
    labs(title = "Disulfide Bonds", x = "Residue Position", y = "Bond Index") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.title.y = element_text(size = 12))
  
  return(ggplotly(p))
}

# Create plots
p1_interactive <- plot_pLDDT(g2p_data)
p2_interactive <- plot_ASA(g2p_data)
p3_interactive <- plot_disulfide_bonds(disulfide_pairs)

# Arrange all plots in a vertical layout
subplot(p1_interactive, p2_interactive, p3_interactive,
        nrows = 3, shareX = TRUE) %>%
  layout(showlegend = TRUE)
